#!/bin/bash

# Set the dynamic linker path to include the directories of Pythia8 and HepMC3 libraries
# export DYLD_LIBRARY_PATH=/opt/homebrew/Cellar/pythia/8.310/lib:/opt/homebrew/Cellar/hepmc3/3.2.7_1/lib:$DYLD_LIBRARY_PATH
#!/usr/bin/env bash
export CC=clang
export CXX=clang++

# (optional) if you use brew-installed HepMC3 and Pythia8:
HPM3_PREFIX=$(brew --prefix hepmc3)
PY8_CXXFLAGS=$(pythia8-config --cxxflags)
PY8_LDFLAGS=$(pythia8-config --ldflags)

clang++ -std=c++17 \
    $PY8_CXXFLAGS \
    -I${HPM3_PREFIX}/include \
    generator.cpp \
    -L${HPM3_PREFIX}/lib -lhepmc3 \
    $PY8_LDFLAGS \
    -o generator
    
# Compile the generator program
# clang++ --std=c++17 \
#     -o generator generator.cpp \
#     -I/opt/homebrew/Cellar/pythia/8.310/include \
#     -I/opt/homebrew/Cellar/hepmc3/3.2.7_1/build/include \
#     -L/opt/homebrew/Cellar/pythia/8.310/lib \
#     -L/opt/homebrew/Cellar/hepmc3/3.2.7_1/build/lib \
#     -lpythia8 -lHepMC3

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful, running the generator..."
else
    echo "Compilation failed"
fi