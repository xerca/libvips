SUBDIRS = overrides

EXTRA_DIST = \
	__init__.py

install-data-hook:
	-if [ ! -f $(pyexecdir)/gi/__init__.py ]; then \
		$(INSTALL_DATA) __init__.py $(pyexecdir)/gi; \
	fi

