/* shrink with a box filter
 *
 * Copyright: 1990, N. Dessipris.
 *
 * Authors: Nicos Dessipris and Kirk Martinez
 * Written on: 29/04/1991
 * Modified on: 2/11/92, 22/2/93 Kirk Martinez - Xres Yres & cleanup 
 incredibly inefficient for box filters as LUTs are used instead of + 
 Needs converting to a smoother filter: eg Gaussian!  KM
 * 15/7/93 JC
 *	- rewritten for partial v2
 *	- ANSIfied
 *	- now shrinks any non-complex type
 *	- no longer cloned from im_convsub()
 *	- could be much better! see km comments above
 * 3/8/93 JC
 *	- rounding bug fixed
 * 11/1/94 JC
 *	- problems with .000001 and round up/down ignored! Try shrink 3738
 *	  pixel image by 9.345000000001
 * 7/10/94 JC
 *	- IM_NEW and IM_ARRAY added
 *	- more typedef
 * 3/7/95 JC
 *	- IM_CODING_LABQ handling added here
 * 20/12/08
 * 	- fall back to im_copy() for 1/1 shrink
 * 2/2/11
 * 	- gtk-doc
 * 10/2/12
 * 	- shrink in chunks to reduce peak memuse for large shrinks
 * 	- simpler
 * 12/6/12
 * 	- redone as a class
 * 	- warn about non-int shrinks
 * 	- some tuning .. tried an int coordinate path, not worthwhile
 */

/*

    This file is part of VIPS.
    
    VIPS is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */

/*

    These files are distributed with VIPS - http://www.vips.ecs.soton.ac.uk

 */

/*
#define VIPS_DEBUG
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/
#include <vips/intl.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vips/vips.h>
#include <vips/debug.h>

#include "resample.h"

typedef struct _VipsShrink {
	VipsResample parent_instance;

	double xshrink;		/* Shrink factors */
	double yshrink;

	/* Size of area of input we average for each output pixel, in pixels.
	 */
	int mw;
	int mh;

	/* Number of pels we average.
	 */
	int np;

	/* We shrink to this, then tilecache to output. 
	 *
	 * We use the tilecache to singlethread requests, then thread again
	 * when we subdivide the input.
	 */
	VipsImage *t;

} VipsShrink;

typedef VipsResampleClass VipsShrinkClass;

G_DEFINE_TYPE( VipsShrink, vips_shrink, VIPS_TYPE_RESAMPLE );

/* Make a sequence value.
 */
static void *
vips_shrink_start( VipsImage *out, void *a, void *b )
{
	VipsImage *in = (VipsImage *) a;

	return( (void *) VIPS_ARRAY( NULL, in->Bands, double ) ); 
}

/* Integer shrink. 
 */
#define ISHRINK( TYPE ) { \
	int *sum = (int *) pels; \
	TYPE *p = (TYPE *) in; \
	TYPE *q = (TYPE *) out; \
	\
	for( j = 0; j < bands; j++ ) \
		sum[j] = 0; \
	\
	for( y1 = 0; y1 < shrink->mh; y1++ ) { \
		for( i = 0, x1 = 0; x1 < shrink->mw; x1++ ) \
			for( j = 0; j < bands; j++, i++ ) \
				sum[j] += p[i]; \
		\
		p += ls; \
	} \
	\
	for( j = 0; j < bands; j++ ) \
		q[j] = (sum[j] + shrink->np / 2) / shrink->np; \
} 

/* Float shrink. 
 */
#define FSHRINK( TYPE ) { \
	double *sum = (double *) pels; \
	TYPE *p = (TYPE *) in; \
	TYPE *q = (TYPE *) out; \
	\
	for( j = 0; j < bands; j++ ) \
		sum[j] = 0.0; \
	\
	for( y1 = 0; y1 < shrink->mh; y1++ ) { \
		for( i = 0, x1 = 0; x1 < shrink->mw; x1++ ) \
			for( j = 0; j < bands; j++, i++ ) \
				sum[j] += p[i]; \
		\
		p += ls; \
	} \
	\
	for( j = 0; j < bands; j++ ) \
		q[j] = sum[j] / shrink->np; \
} 

/* Have one of these for each call to vips_shrink_gen().
 */
typedef struct _GenArgs {
	VipsShrink *shrink;
	VipsRegion *or;
} GenArgs;

static int
vips_shrink_gen2( VipsRegion *ir, 
	void *seq, void *a, void *b, gboolean *stop )
{
	VipsPel *pels = (VipsPel *) seq;
	VipsImage *im = (VipsImage *) a;
	GenArgs *args = (GenArgs *) b;
	VipsShrink *shrink = args->shrink;
	VipsRegion *or = args->or;
	VipsRect *r = &ir->valid;
	const int bands = im->Bands;
	const int sizeof_pixel = VIPS_IMAGE_SIZEOF_PEL( im );
	const int ls = VIPS_REGION_LSKIP( ir ) / 
		VIPS_IMAGE_SIZEOF_ELEMENT( im );

	int x, y, i;
	int x1, y1, j;

	int left;
	int top;
	int width;
	int height;

	/*
	VIPS_DEBUG_MSG( "vips_shrink_gen2: %p left = %d, top = %d, "
		"width = %d, height = %d\n", 
		seq,
		r->left, r->top, r->width, r->height ); 
	 */

	/* Corresponding output rect.
	 */
	left = r->left / shrink->mw;
	top = r->top / shrink->mh;
	width = r->width / shrink->mw;
	height = r->height / shrink->mh;

	for( y = 0; y < height; y++ ) { 
		VipsPel *out = VIPS_REGION_ADDR( or, left, top + y ); 

		for( x = 0; x < width; x++ ) { 
			int ix = (left + x) * shrink->xshrink; 
			int iy = (top + y) * shrink->yshrink; 

			VipsPel *in = VIPS_REGION_ADDR( ir, ix, iy ); 

			switch( im->BandFmt ) {
			case VIPS_FORMAT_UCHAR: 	
				ISHRINK( unsigned char ); break;
			case VIPS_FORMAT_CHAR: 	
				ISHRINK( char ); break; 
			case VIPS_FORMAT_USHORT: 
				ISHRINK( unsigned short ); break;
			case VIPS_FORMAT_SHORT: 	
				ISHRINK( short ); break; 
			case VIPS_FORMAT_UINT: 	
				ISHRINK( unsigned int ); break; 
			case VIPS_FORMAT_INT: 	
				ISHRINK( int );  break; 
			case VIPS_FORMAT_FLOAT: 	
				FSHRINK( float ); break; 
			case VIPS_FORMAT_DOUBLE:	
				FSHRINK( double ); break;

			default:
				g_assert( 0 ); 
			}

			out += sizeof_pixel;
		}
	}

	return( 0 );
}

/* Free a sequence value.
 */
static int
vips_shrink_stop( void *vseq, void *a, void *b )
{
	g_free( vseq ); 

	return( 0 );
}

static int
vips_shrink_gen( VipsRegion *or, void *vseq, void *a, void *b, gboolean *stop )
{
	VipsImage *in = (VipsImage *) a;
	VipsShrink *shrink = (VipsShrink *) b;

	/* How do we step through the input image? 
	 *
	 * We want to do it in chunks which are multiples of mw/mh, so that we
	 * can average each chunk to make a whole number of outputs. We want
	 * to stay around VIPS__TILE_WIDTH/VIPS__TILE_HEIGHT pixels for each
	 * chunk.
	 */
	int tile_width = shrink->mw > VIPS__TILE_WIDTH ?
		shrink->mw : shrink->mw * (VIPS__TILE_WIDTH / shrink->mw);
	int tile_height = shrink->mh > VIPS__TILE_HEIGHT ?
		shrink->mh : shrink->mh * (VIPS__TILE_HEIGHT / shrink->mh);

	GenArgs args;
	VipsRect area;

	VIPS_DEBUG_MSG( "vips_shrink_gen: thread = %p, left = %d, top = %d, "
		"width = %d, height = %d\n", 
		g_thread_self(),
		or->valid.left, or->valid.top, 
		or->valid.width, or->valid.height ); 

	args.shrink = shrink;
	args.or = or;

	/* We need this area in our source image.
	 */
	area.left = or->valid.left * shrink->xshrink;
	area.top = or->valid.top * shrink->yshrink;
	area.width = ceil( or->valid.width * shrink->xshrink );
	area.height = ceil( or->valid.height * shrink->yshrink );

	/*
	VIPS_DEBUG_MSG( "vips_shrink_gen: scanning left = %d, top = %d, "
		"width = %d, height = %d\n", 
		area.left, area.top, area.width, area.height ); 
	VIPS_DEBUG_MSG( "vips_shrink_gen: tile_width = %d, tile_height = %d\n",
		tile_width, tile_height ); 
	 */

	/* Scan that chunk of our source image with threads.
	 */
	if( vips_sink_area( in, 
		tile_width, tile_height, 
		&area,
		vips_shrink_start, vips_shrink_gen2, vips_shrink_stop, 
		in, &args ) )
		return( -1 ); 

	return( 0 );
}

static int
vips_shrink_build( VipsObject *object )
{
	VipsResample *resample = VIPS_RESAMPLE( object );
	VipsShrink *shrink = (VipsShrink *) object;

	VipsImage *x;

	if( VIPS_OBJECT_CLASS( vips_shrink_parent_class )->build( object ) )
		return( -1 );

	shrink->mw = ceil( shrink->xshrink );
	shrink->mh = ceil( shrink->yshrink );
	shrink->np = shrink->mw * shrink->mh;

	if( im_check_noncomplex( "VipsShrink", resample->in ) )
		return( -1 );

	if( shrink->xshrink < 1.0 || 
		shrink->yshrink < 1.0 ) {
		vips_error( "VipsShrink", 
			"%s", _( "shrink factors should be >= 1" ) );
		return( -1 );
	}

	if( (int) shrink->xshrink != shrink->xshrink || 
		(int) shrink->yshrink != shrink->yshrink ) 
		vips_warn( "VipsShrink", 
			"%s", _( "not integer shrink factors, "
				"expect poor results" ) ); 

	if( shrink->xshrink == 1.0 &&
		shrink->yshrink == 1.0 )
		return( vips_image_write( resample->in, resample->out ) );

	/* Shrink to this, tilecache to output.
	 */
	shrink->t = vips_image_new();
	vips_object_local( shrink, shrink->t );

	if( vips_image_copy_fields( shrink->t, resample->in ) )
		return( -1 );

	/* THINSTRIP will work, FATSTRIP will break seq mode. If you combine
	 * shrink with conv you'll need to use a line cache to maintain
	 * sequentiality.
	 */
	vips_demand_hint( shrink->t, 
		VIPS_DEMAND_STYLE_THINSTRIP, resample->in, NULL );

	/* Size output. Note: we round the output width down!
	 */
	shrink->t->Xsize = resample->in->Xsize / shrink->xshrink;
	shrink->t->Ysize = resample->in->Ysize / shrink->yshrink;
	shrink->t->Xres = resample->in->Xres / shrink->xshrink;
	shrink->t->Yres = resample->in->Yres / shrink->yshrink;
	if( shrink->t->Xsize <= 0 || 
		shrink->t->Ysize <= 0 ) {
		vips_error( "VipsShrink", 
			"%s", _( "image has shrunk to nothing" ) );
		return( -1 );
	}

	if( vips_image_generate( shrink->t,
		NULL, vips_shrink_gen, NULL, 
		resample->in, shrink ) )
		return( -1 );

	if( vips_tilecache( shrink->t, &x, 
		"tile_width", VIPS__TILE_WIDTH,
		"tile_height", 1,
		"max_tiles", 2 * VIPS__TILE_HEIGHT * 
			((shrink->t->Xsize / VIPS__TILE_WIDTH) + 1),
		NULL ) )
		return( -1 ); 

	g_object_set( shrink, "out", x, NULL );

	return( 0 );
}

static void
vips_shrink_class_init( VipsShrinkClass *class )
{
	GObjectClass *gobject_class = G_OBJECT_CLASS( class );
	VipsObjectClass *vobject_class = VIPS_OBJECT_CLASS( class );
	VipsOperationClass *operation_class = VIPS_OPERATION_CLASS( class );

	VIPS_DEBUG_MSG( "vips_shrink_class_init\n" );

	gobject_class->set_property = vips_object_set_property;
	gobject_class->get_property = vips_object_get_property;

	vobject_class->nickname = "shrink";
	vobject_class->description = _( "shrink an image" );
	vobject_class->build = vips_shrink_build;

	operation_class->flags = VIPS_OPERATION_SEQUENTIAL;

	VIPS_ARG_DOUBLE( class, "xshrink", 8, 
		_( "Xshrink" ), 
		_( "Horizontal shrink factor" ),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET( VipsShrink, xshrink ),
		1.0, 1000000, 1 );

	VIPS_ARG_DOUBLE( class, "yshrink", 9, 
		_( "Yshrink" ), 
		_( "Vertical shrink factor" ),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET( VipsShrink, yshrink ),
		1.0, 1000000, 1 );

}

static void
vips_shrink_init( VipsShrink *shrink )
{
}

/**
 * vips_shrink:
 * @in: input image
 * @out: output image
 * @xshrink: horizontal shrink
 * @yshrink: vertical shrink
 *
 * Shrink @in by a pair of factors with a simple box filter. 
 *
 * You will get aliasing for non-integer shrinks. In this case, shrink with
 * this function to the nearest integer size above the target shrink, then
 * downsample to the exact size with im_affinei() and your choice of
 * interpolator.
 *
 * See also: im_affinei().
 *
 * Returns: 0 on success, -1 on error
 */
int
vips_shrink( VipsImage *in, VipsImage **out, 
	double xshrink, double yshrink, ... )
{
	va_list ap;
	int result;

	va_start( ap, yshrink );
	result = vips_call_split( "shrink", ap, in, out, xshrink, yshrink );
	va_end( ap );

	return( result );
}
