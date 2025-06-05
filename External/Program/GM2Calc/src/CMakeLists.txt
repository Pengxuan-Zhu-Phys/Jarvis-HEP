# GM2Calc library
add_library(gm2calc
  gm2_dilog.cpp
  gm2_error_c.cpp
  gm2_ffunctions.cpp
  gm2_mf.cpp
  gm2_numerics.cpp
  gm2_slha_io.cpp
  MSSMNoFV/gm2_1loop_c.cpp
  MSSMNoFV/gm2_1loop.cpp
  MSSMNoFV/gm2_2loop_c.cpp
  MSSMNoFV/gm2_2loop.cpp
  MSSMNoFV/gm2_uncertainty_c.cpp
  MSSMNoFV/gm2_uncertainty.cpp
  MSSMNoFV/MSSMNoFV_onshell_c.cpp
  MSSMNoFV/MSSMNoFV_onshell.cpp
  MSSMNoFV/MSSMNoFV_onshell_mass_eigenstates.cpp
  MSSMNoFV/MSSMNoFV_onshell_physical.cpp
  MSSMNoFV/MSSMNoFV_onshell_problems.cpp
  MSSMNoFV/MSSMNoFV_onshell_soft_parameters.cpp
  MSSMNoFV/MSSMNoFV_onshell_susy_parameters.cpp
  THDM/gm2_1loop.cpp
  THDM/gm2_1loop_c.cpp
  THDM/gm2_1loop_H.cpp
  THDM/gm2_2loop.cpp
  THDM/gm2_2loop_c.cpp
  THDM/gm2_2loop_B.cpp
  THDM/gm2_2loop_F.cpp
  THDM/gm2_uncertainty.cpp
  THDM/gm2_uncertainty_c.cpp
  THDM/THDM.cpp
  THDM/THDM_c.cpp
  THDM/THDM_mass_eigenstates.cpp
  THDM/THDM_problems.cpp
  THDM/THDM_parameters.cpp
  SM/SM.cpp
  SM/SM_c.cpp
  )
target_include_directories(gm2calc
  PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_link_libraries(gm2calc
  PUBLIC
    Eigen3::Eigen
  PRIVATE
    Boost::boost
)
add_library(GM2Calc::GM2Calc ALIAS gm2calc)

install(
  TARGETS gm2calc
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# GM2Calc main executable
add_executable(gm2calc.x gm2calc.cpp)
target_link_libraries(gm2calc.x PRIVATE GM2Calc::GM2Calc)
target_include_directories(gm2calc.x
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${Boost_INCLUDE_DIRS}"
)

install(
  TARGETS gm2calc.x
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# MathML executable
if(Mathematica_MathLink_FOUND)
  Mathematica_MathLink_ADD_EXECUTABLE(
    gm2calc.mx
    gm2calc.tm)
  target_link_libraries(gm2calc.mx PRIVATE GM2Calc::GM2Calc)
  target_include_directories(gm2calc.mx
    PRIVATE
      "${Mathematica_MathLink_INCLUDE_DIR}"
      "${CMAKE_CURRENT_SOURCE_DIR}"
      "${Boost_INCLUDE_DIRS}"
  )
  install(
    TARGETS gm2calc.mx
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
