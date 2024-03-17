#!/bin/csh
if( ! $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH /Users/buding/Workshop/Jarvis/External/Library/HepMC/lib
else
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/Users/buding/Workshop/Jarvis/External/Library/HepMC/lib
endif
setenv PYTHIA8DATA ${PYTHIA8_HOME}/xmldoc
