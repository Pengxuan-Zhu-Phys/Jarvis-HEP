CC = gcc
FC = gfortran
#CXX =
#CFLAGS = 
#FFLAGS = 
#CXXFLAGS = 


.PHONY: all clean


all: libSLHAplus.a

OBJ =  alpha_s.o aux.o ferror.o  jacobi.o fixArg.o chDiag.o SLHAreader.o 



libSLHAplus.a:libSLHAplus.a($(OBJ))
	rm -f *.o

clean: 
	rm -f libSLHAplus.a *.o
