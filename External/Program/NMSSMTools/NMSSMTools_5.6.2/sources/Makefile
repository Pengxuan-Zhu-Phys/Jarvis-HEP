.PHONY: clean

ifneq ($(MAKECMDGOALS),clean)
include micromegas/CalcHEP_src/FlagsForMake
endif

all:	sources.a 

SRC = $(wildcard *.f) 
OBJ = $(SRC:.f=.o) relden.o

sources.a:sources.a($(OBJ)) 
	$(RANLIB) sources.a 
	rm -f *.o

clean:
	rm -f *.a  *.o
