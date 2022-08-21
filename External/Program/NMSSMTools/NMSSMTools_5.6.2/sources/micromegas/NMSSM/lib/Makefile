.PHONY: clean all

NMSSMTools = $(CURDIR)/../../Packages/NMSSMTools_5.5.2
CALCHEP =    ../../CalcHEP_src
all:  aLib.a $(NMSSMTools)/sources/sources.a  $(NMSSMTools)/main/nmhdecay  $(NMSSMTools)/main/nmspec

ifneq ($(MAKECMDGOALS),clean)
include $(CALCHEP)/FlagsForMake
endif

aLibOBJ = $(notdir $(subst .c,.o,$(wildcard *.c))) fort.o


aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)

aLib.a(calchep_func.o):lpath.h

lpath.h:
	echo \#define LPATH \"`pwd`\"  >lpath.h
	echo \#define NMSSMTOOLS \"$(NMSSMTools)\" >> lpath.h

$(NMSSMTools)/sources/sources.a:
	make -C $(NMSSMTools) init

$(NMSSMTools)/main/nmhdecay: $(NMSSMTools)/sources/sources.a
	make -C $(NMSSMTools)/main  nmhdecay

$(NMSSMTools)/main/nmspec: $(NMSSMTools)/sources/sources.a
	make -C $(NMSSMTools)/main nmspec



clean: 
	rm -f lpath.h *.o *.a  
	make -C nngg-nm4b clean