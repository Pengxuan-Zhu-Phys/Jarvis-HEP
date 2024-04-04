
set( example_code 
		  example_MyHerwig.cc
		  example_MyPythia.cc
		  example_MyPythiaOnlyToHepMC.cc
		  example_PythiaStreamIO.cc
		  initpydata.f
		  initPythia.cc
		  PythiaHelper.h
		  testHerwigCopies.cc
		  testPythiaCopies.cc
		    )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/GNUmakefile.example.in 
                ${CMAKE_CURRENT_BINARY_DIR}/GNUmakefile  @ONLY )

INSTALL (FILES ${example_code}  
               ${CMAKE_CURRENT_BINARY_DIR}/GNUmakefile
         DESTINATION share/HepMC/examples/fio )
