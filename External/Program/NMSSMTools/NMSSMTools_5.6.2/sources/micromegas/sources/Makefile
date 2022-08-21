  
.PHONY: all clean

ifneq ($(MAKECMDGOALS),clean)
include ../CalcHEP_src/FlagsForMake
endif

CFLAGS+= -I../include
all:../lib/micromegas.a ../lib/dummy.a ../lib/maxGap.so

special=improveCS.c

OBJ = $(notdir $(subst .c,.o,$(filter-out $(special),$(wildcard *.c))))


#ifneq ($(FC),)
#   override OBJ += faux.o
#endif
../lib/micromegas.a(cs22.o):data/et_tab.inc

../lib/micromegas.a:../lib/micromegas.a($(OBJ)) 
	$(RANLIB) ../lib/micromegas.a 
	rm -f *.o

../lib/dummy.a:../lib/dummy.a(improveCS.o)
	rm -f improveCS.o
	
../lib/maxGap.so:
ifneq ($(FC),)
	./maxGap/makePath $(CURDIR)/maxGap >maxGap/pathData.f
	$(FC) $(FFLAGS)  -shared  -o ../lib/maxGap.so -ffixed-line-length-none -fPIC maxGap/*.f 	
else
	$(CC) $(CFLAGS)  -shared  -o ../lib/maxGap.so -fPIC maxGap/maxGapDummy.c 
endif

clean:
	rm -f ../lib/*.a ../lib/*.so  *.o 
