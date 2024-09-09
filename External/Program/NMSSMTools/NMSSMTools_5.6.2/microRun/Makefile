

.PHONY: all libs clean specialClean

.PHONY: all libs clean specialClean

# read compiler flags
ifneq ($(MAKECMDGOALS),clean)
AllFlags = ../sources/micromegas/CalcHEP_src/FlagsForMake
ifeq (,$(wildcard $(AllFlags) )) 
$(error File $(AllFlags) is absent. Presumably you forgot to compile main code)
endif 
include ../sources/micromegas/CalcHEP_src/FlagsForMake
endif



cLib = $(CALCHEP)/lib
# files to compile
SSS = $(wildcard lib/*.a) ../sources/micromegas/lib/micromegas.a  $(cLib)/dynamic_me.a ../sources/micromegas/lib/micromegas.a \
 work/work_aux.a  $(wildcard lib/*.a)   $(cLib)/sqme_aux.$(SO) $(cLib)/libSLHAplus.a   \
 $(cLib)/num_c.a   $(cLib)/serv.a $(cLib)/ntools.a   ../sources/micromegas/lib/maxGap.so  ../sources/micromegas/lib/dummy.a   $(LX11)


ifneq ($(LHAPDFPATH),)
  SSS += -L$(LHAPDFPATH)  -lLHAPDF $(cLib)/dummy.a
  DLSET = export LD_RUN_PATH=$(LHAPDFPATH);
else 
  SSS += $(cLib)/dummy.a 
  DLSET= 
endif   


ifeq (,$(main)) 
  $(warning Main program is not specified. Use gmake main='<code of main program>')
  $(warning Default  main=main.c)
  main = main.c
endif


all: $(main) libs calchep/bin
  ifeq (.cpp,$(suffix $(main)))
	$(DLSET) $(CXX) $(CXXFLAGS) -o $(main:.cpp=) $(main)  $(SSS)  $(lDL)  -lm -lpthread
  else  
    ifeq (.c,$(suffix $(main)))
	$(DLSET) $(CC)  $(CFLAGS) -o $(main:.c=) $(main)  $(SSS) $(lDL)  -lm -lpthread
    else  
      ifeq (.F,$(suffix $(main)))
        ifeq ($(FC), )
	$(error  Fortran compiler was not detected.)
        else 
	$(DLSET) $(FC) $(FFLAGS) -ffixed-line-length-none -o $(main:.F=) $(main) $(cLib)/faux.o  $(SSS)   $(lDL)  -lm -lpthread
        endif
      else
	$(error Only .c, .cpp, and .F  main programs  are supported).
      endif
    endif
  endif  



libs:
	$(MAKE) -C work
	$(MAKE) -C lib
	$(MAKE) -C ../sources
calchep/bin:
	ln -s  `pwd`/../sources/micromegas/CalcHEP_src/bin  `pwd`/calchep/bin

clean:  specialClean
	../sources/micromegas/sources/cleanexe
	rm -rf *.dSYM
	$(MAKE) -C lib  clean
	$(MAKE) -C work clean
	@-unlink calchep/bin
	rm -rf calchep/tmp/* calchep/results/* 
	rm -f HB.* HS.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat
	rm -f Lilith_*   particles.py*
	rm -f  smodels.* summary.*



specialClean:
	rm -f  lanhep/*.mdl lanhep/masses.chk
	rm -f  omega decay spectr inp nngg.out
	rm -f  lib/nngg-nm4b/lGamma.exe
