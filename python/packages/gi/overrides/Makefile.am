vips_overridesdir = $(pyexecdir)/gi/overrides
vips_overrides_PYTHON = Vips.py

EXTRA_DIST = \
	__init__.py

install-data-hook:
	if [ ! -f $(pyexecdir)/gi/overrides/__init__.py ]; then \
		$(INSTALL_DATA) __init__.py $(pyexecdir)/gi/overrides; \
	fi

