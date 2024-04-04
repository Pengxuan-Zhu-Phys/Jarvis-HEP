.PHONY: clean

LOOPTOOLS = ../../../Packages/LoopTools-2.12-v0/

lGamma.exe: $(LOOPTOOLS)/FlagsForMake $(LOOPTOOLS)/build/libooptools.a
	$(MAKE)  -f makeCompile
	
 $(LOOPTOOLS)/FlagsForMake:
	cd $(LOOPTOOLS); ./configure
 	
 $(LOOPTOOLS)/build/libooptools.a: $(LOOPTOOLS)/FlagsForMake
	 $(MAKE)  -C $(LOOPTOOLS)

clean:
	rm -f nngg.a util.a lGamma.exe *.o */*.o 	