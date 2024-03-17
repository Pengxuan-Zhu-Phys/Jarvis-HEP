# This makefile is used to copy files during the "make install" step

INSTALLDIR = $(pkgdatadir)/examples/pythia8

# files to distribute
EXTRA_DIST = \
    main31.cc \
    main32.cc \
    main32.cmnd \
    README

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(INSTALLDIR)
	$(INSTALL_DATA) GNUmakefile.example $(DESTDIR)$(INSTALLDIR)/GNUmakefile
	$(INSTALL_DATA) config.csh $(DESTDIR)$(INSTALLDIR)/config.csh
	$(INSTALL_DATA) config.sh $(DESTDIR)$(INSTALLDIR)/config.sh
	for file in $(EXTRA_DIST); do \
	  $(INSTALL_DATA) $(srcdir)/$$file $(DESTDIR)$(INSTALLDIR)/$$file; done
