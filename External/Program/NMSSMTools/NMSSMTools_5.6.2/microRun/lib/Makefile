.PHONY: clean all

Tools = ../../sources/micromegas/Packages/LoopTools-2.12-v0

ifneq ($(MAKECMDGOALS),clean)
include ../../sources/micromegas/CalcHEP_src/FlagsForMake
endif

aLibOBJ = $(notdir $(subst .c,.o,$(wildcard *.c))) fort.o


all:  aLib.a $(Tools)/build/libooptools.a  nngg-nm4b/lGamma.exe

$(Tools)/build/libooptools.a:
	cd $(Tools); ./configure
	make -C $(Tools)
	
nngg-nm4b/lGamma.exe: $(Tools)/build/libooptools.a
	 make -C nngg-nm4b

aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)

aLib.a(calchep_func.o):lpath.h

clean: 
	rm -f lpath.h *.o *.a  
	make -C nngg-nm4b  clean 
	make -C $(Tools) clean
