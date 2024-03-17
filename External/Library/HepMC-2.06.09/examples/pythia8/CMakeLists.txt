
set( example_code 
		  main31.cc
		  main32.cc
		  main32.cmnd
		  README
		  )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/GNUmakefile.example.in 
                ${CMAKE_CURRENT_BINARY_DIR}/GNUmakefile  @ONLY )
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/config.csh.in
                ${CMAKE_CURRENT_BINARY_DIR}/config.csh  @ONLY )
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/config.sh.in
                ${CMAKE_CURRENT_BINARY_DIR}/config.sh  @ONLY )

INSTALL (FILES ${example_code}  
               ${CMAKE_CURRENT_BINARY_DIR}/GNUmakefile
               ${CMAKE_CURRENT_BINARY_DIR}/config.csh
               ${CMAKE_CURRENT_BINARY_DIR}/config.sh
         DESTINATION share/HepMC/examples/pythia8 )
