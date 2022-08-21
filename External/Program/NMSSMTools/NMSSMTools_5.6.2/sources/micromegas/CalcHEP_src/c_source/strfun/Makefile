.PHONY: all

include ../../FlagsForMake

Idir = .. chep_crt plot service2 getmem dynamicME SLHAplus ntools num

override CFLAGS += -Iinclude $(patsubst %,-I../%/include ,$(Idir))

special = pdf_dummy.c lhapdf2pdt.c

OBJ = $(notdir $(subst .c,.o,$(filter-out $(special),$(wildcard *.c ))))

#VPATH = strfun

lib = $(CALCHEP)/lib
bin = ../../bin

all: $(lib)/num_c.a $(lib)/dummy.a $(bin)/lhapdf2pdt

$(lib)/num_c.a:$(lib)/num_c.a($(OBJ))
	$(RANLIB) $(lib)/num_c.a
	@rm -f $(OBJ)

$(lib)/dummy.a:$(lib)/dummy.a(pdf_dummy.o)
	$(RANLIB) $@
	@rm -f pdf_dummy.o 

$(bin)/lhapdf2pdt: lhapdf2pdt.c 
	$(CC) $(CFLAGS)  -o $(bin)/lhapdf2pdt lhapdf2pdt.c $(lib)/ntools.a -lm

