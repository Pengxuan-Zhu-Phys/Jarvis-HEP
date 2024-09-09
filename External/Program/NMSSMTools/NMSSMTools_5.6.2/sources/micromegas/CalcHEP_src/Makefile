
.PHONY: all  COMPILE clean

ifeq ($(MAKECMDGOALS),clean)
clean :
	if(test -r FlagsForSh) then echo FlagsForSh - compiler flags; rm -i  FlagsForSh; fi
	./sbin/setPath " "
	@rm -f  include/rootDir.h
	@rm -f  FlagsForMake  so_locations CMessage
	@cd lib;  rm -rf sqme_aux.so.dSYM;  rm -f *.*
	@cd bin;  rm -rf *.dSYM rm -f Int calc events2tab s_calchep show_distr sum_distr plot_view make_VandP event_info event_mixer lhe2tab nt_maker showHelp event2lhe *.exe lhapdf2pdt
	@rm -rf sbin/make-j sbin/makeVrtLib sbin/*.dSYM   
	@rm -f c_source/*/*.o c_source/*/so_location
	@rm -fr work/results/* work/batch_results/*  work/tmp/*  work/Processes work/Events work/html
	@rm -f c_source/Root/ch_dict.cc c_source/Root/ch_dict.h c_source/Root/*.o
	@-unlink work/bin
	@chmod 644 mkWORKdir
	cp calchep.ini work
endif

all:FlagsForMake COMPILE


flags:
	./getFlags

FlagsForMake: flags

COMPILE:FlagsForMake
	./sbin/setPath $(CURDIR)
	@if( test ! -d work/bin) then ln -s  `pwd`/bin  `pwd`/work/bin; fi
	chmod 755 mkWORKdir
	$(MAKE) -C c_source
	#--------------------------------------------------------
	# CalcHEP has compiled successfuly and can be started.
	# The manual can be found on the CalcHEP website:
	#      http://theory.sinp.msu.ru/~pukhov/calchep.html
	# The next step is typically to run 
	#      ./mkWORKdir  <new_dir>
	# where <new_dir> is the new directory where you will do
	# your calculations.  After creating this directory, you
	# should cd into it and run calchep or calchep_batch.
	# Please see the manual for further details.
	#---------------------------------------------------------"
	
	@if(test -z "`grep lX11 FlagsForMake`") then cat .X11; fi
