# This makefile is used to copy files during the "make install" step

INSTALLDIR = $(pkgdatadir)/examples/fio

# files to distribute
EXTRA_DIST = \
    example_MyHerwig.cc \
    example_MyPythia.cc \
    example_MyPythiaOnlyToHepMC.cc \
    example_PythiaStreamIO.cc \
    testHerwigCopies.cc	\
    testPythiaCopies.cc	\
    initPythia.cc \
    initpydata.f \
    PythiaHelper.h

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(INSTALLDIR)
	$(INSTALL_DATA) GNUmakefile.example $(DESTDIR)$(INSTALLDIR)/GNUmakefile
	for file in $(EXTRA_DIST); do \
	  $(INSTALL_DATA) $(srcdir)/$$file $(DESTDIR)$(INSTALLDIR)/$$file; done
