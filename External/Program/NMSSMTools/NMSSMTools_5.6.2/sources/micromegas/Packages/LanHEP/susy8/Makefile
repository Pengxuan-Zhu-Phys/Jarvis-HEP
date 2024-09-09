CFLAG = -O2
FORTLIB = -lgfortran
SHLIBS = sort4.so sugra.so nmssm.so fhf3.so neutd.so cpviol.so

.PHONY: all clean distclean

all: $(SHLIBS)

sort4.so: sort4.c
	gcc -shared $(CFLAG) -o sort4.so sort4.c

sugra.so: sugra.c
	gcc -shared $(CFLAG) -o sugra.so sugra.c

nmssm.so: nmssm.c
	gcc -shared $(CFLAG) -o nmssm.so nmssm.c

fhf3.so: fhf3.c
	gcc -shared $(CFLAG) -o fhf3.so fhf3.c fhf/*.f $(FORTLIB)
	
neutd.so: neutc.c neutd.F
	gcc -shared $(CFLAG) -o neutd.so neutc.c neutd.F $(FORTLIB)

cpviol.so: cpvc.c cpvf.F
	gcc -shared $(CFLAG) -o cpviol.so cpvc.c cpvf.F $(FORTLIB)

clean:
	rm -f *.in *.out slhaInput slhaInput.sav\
		slhaOutput slhaOutput.sav

distclean: clean
	rm -f  $(SHLIBS)
