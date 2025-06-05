cmake_minimum_required(VERSION 3.12)

project(GM2Calc
  DESCRIPTION "Calculation of BSM contributions to the anomalous magnetic moment of the muon"
  LANGUAGES C CXX)

option(ENABLE_EXAMPLES "Enable examples" ON)
option(ENABLE_MATHEMATICA "Enable Mathematica interface" ON)
option(ENABLE_PYTHON "Enable Python interface" ON)
option(ENABLE_TESTS "Enable tests" ON)

include(CTest)
include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED YES)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${PROJECT_SOURCE_DIR}/cmake/Mathematica)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# determine GM2Calc version
file(
  STRINGS "${PROJECT_SOURCE_DIR}/include/gm2calc/gm2_version.h" _gm2calc_version_line
  REGEX "^#define[ \t]+GM2CALC_VERSION[ \t]+\"([0-9.]+)\"$")
string(REGEX REPLACE "^#define[ \t]+GM2CALC_VERSION[ \t]+\"([0-9.]+)\"$" "\\1" GM2Calc_VERSION "${_gm2calc_version_line}")

find_package(Boost 1.37.0 REQUIRED)
find_package(Eigen3 3.1 REQUIRED NO_MODULE)

# conan's Eigen3 package defines Eigen3_INCLUDE_DIR instead of
# EIGEN3_INCLUDE_DIR
if(NOT DEFINED EIGEN3_INCLUDE_DIR)
  if(DEFINED Eigen3_INCLUDE_DIR)
    set(EIGEN3_INCLUDE_DIR "${Eigen3_INCLUDE_DIR}")
  endif()
endif()

if(ENABLE_MATHEMATICA)
  find_package(Mathematica COMPONENTS MathLink)
endif()

if(ENABLE_PYTHON)
  # Needed a shared library for cppyy package used for the python interface
  if(${BUILD_SHARED_LIBS})
    find_package(Python3 COMPONENTS Interpreter Development)
    set(GM2CALC_VENV "${PROJECT_BINARY_DIR}/gm2calc-venv")
    execute_process(COMMAND "${Python3_EXECUTABLE}" -m venv "${GM2CALC_VENV}" OUTPUT_QUIET)
    set(ENV{VIRTUAL_ENV} "${GM2CALC_VENV}")
    set(Python3_FIND_VIRTUALENV FIRST)
    unset(Python3_EXECUTABLE)
    find_package(Python3 COMPONENTS Interpreter Development)
    # python executable in the virtual environment
    set(PYTHON_EXE "${GM2CALC_VENV}/bin/python")
  else()
    message("-- Python interface requires a shared library to be built.  Use the option -DBUILD_SHARED_LIBS=ON to enable.")
  endif()
endif()

if(Python3_FOUND)
  execute_process(
    COMMAND "${PYTHON_EXE}" -m pip install cppyy
    OUTPUT_QUIET)
  execute_process(
    COMMAND "${PYTHON_EXE}" -c "import cppyy"
    RESULT_VARIABLE PYTHON3_HAS_cppyy
    OUTPUT_QUIET)
  if(NOT ${PYTHON3_HAS_cppyy} EQUAL 0)
    message("-- Python3 does not have cppyy.  cppyy module is required for python interface.")
  else()
    set(Python_INTERFACE 1)
  endif()
endif()

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(doc)

if(ENABLE_TESTS AND EXISTS "${PROJECT_SOURCE_DIR}/test")
  enable_testing()
  add_subdirectory(test)
endif()

# install public headers
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/include/gm2calc"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

# generate gm2calc.pc file for pkg-config
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/gm2calc.pc.in"
  "${PROJECT_BINARY_DIR}/gm2calc.pc"
  @ONLY)

# install gm2calc.pc file
install(
  FILES "${PROJECT_BINARY_DIR}/gm2calc.pc"
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig")

if(${Python_INTERFACE})
  message("-- Copying python examples")
  configure_file(${PROJECT_SOURCE_DIR}/src/gm2_python_interface.py.in 
                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gm2_python_interface.py 
                 @ONLY)
  configure_file(${PROJECT_SOURCE_DIR}/examples/example_thdm.py.in 
                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/example_thdm.py 
                 @ONLY)
  configure_file(${PROJECT_SOURCE_DIR}/examples/example_slha.py.in
                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/example_slha.py 
                 @ONLY)
  configure_file(${PROJECT_SOURCE_DIR}/examples/example_gm2calc.py.in 
                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/example_gm2calc.py 
                 @ONLY)
  configure_file(${PROJECT_SOURCE_DIR}/examples/example_gm2scan.py.in
                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/example_gm2scan.py 
                 @ONLY)
endif()

