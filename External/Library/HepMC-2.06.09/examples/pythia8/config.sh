#!/bin/sh
if [ ! $?LD_LIBRARY_PATH ]; then
  export LD_LIBRARY_PATH=/Users/buding/Workshop/Jarvis/External/Library/HepMC/lib
fi
if [ $?LD_LIBRARY_PATH ]; then
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/Users/buding/Workshop/Jarvis/External/Library/HepMC/lib
fi
export PYTHIA8DATA=${PYTHIA8_HOME}/xmldoc
