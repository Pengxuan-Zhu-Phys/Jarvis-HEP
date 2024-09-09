
.PHONY: all clean flags

Tools = Packages/LoopTools-2.12-v0

all:include/microPath.h $(Tools)/build/libooptools.a
	$(MAKE) -C CalcHEP_src  MICROMEGAS=MICROMEGAS
	$(MAKE) -C sources
	
$(Tools)/build/libooptools.a:
	cd $(Tools);  ./configure
	make -C $(Tools)
        
include/microPath.h:
	echo \#define micrO \"$(CURDIR)\"  > include/microPath.h
        
   
clean:  
	rm -f include/microPath.h
	./clean
	make -C $(Tools)  clean
	rm -rf */*.dSYM 

flags: 
	$(MAKE) -C CalcHEP_src flags