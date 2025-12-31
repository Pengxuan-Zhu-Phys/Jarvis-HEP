# Distributed under the GPLv3 License.
# Author: Alexander Voigt
#
# FindGM2Calc
# -----------
#
# Finds the GM2Calc library.
# [https://github.com/GM2Calc/GM2Calc]
#
# This module reads the following variables:
#
# GM2Calc_LIBRARY       - GM2Calc library directory
# GM2Calc_INCLUDE_DIR   - GM2Calc include directory
#
# This module defines the following variables:
#
# GM2Calc_FOUND         - set if GM2Calc has been found
# GM2Calc_VERSION       - GM2Calc version
# GM2Calc_INCLUDE_DIRS  - GM2Calc include directory
# GM2Calc_LIBRARIES     - GM2Calc library
#
# and defines the following imported targets:
#
# GM2Calc::GM2Calc

# search gm2calc/gm2_version.h first in ${GM2Calc_INCLUDE_DIR}
find_path(GM2Calc_INCLUDE_DIRS
  NAMES gm2calc/gm2_version.h
  PATHS
    ${GM2Calc_INCLUDE_DIR}
  PATH_SUFFIXES
    include
  NO_DEFAULT_PATH
  NO_CMAKE_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  NO_CMAKE_FIND_ROOT_PATH
)

find_path(GM2Calc_INCLUDE_DIRS
  NAMES gm2calc/gm2_version.h
)

# search GM2Calc library first in ${GM2Calc_LIBRARY}
find_library(GM2Calc_LIBRARIES
  NAMES gm2calc
  PATHS
    ${GM2Calc_LIBRARY}
  PATH_SUFFIXES
    build
  NO_DEFAULT_PATH
  NO_CMAKE_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  NO_CMAKE_FIND_ROOT_PATH
)

find_library(GM2Calc_LIBRARIES
  NAMES gm2calc
)

# find version
if(GM2Calc_INCLUDE_DIRS)
  file(READ "${GM2Calc_INCLUDE_DIRS}/gm2calc/gm2_version.h" _gm2calc_version_header)

  string(REGEX MATCH "define[ \t]+GM2CALC_VERSION_MAJOR[ \t]+([0-9]+)" _gm2calc_version_major_match "${_gm2calc_version_header}")
  set(GM2Calc_VERSION_MAJOR "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+GM2CALC_VERSION_MINOR[ \t]+([0-9]+)" _gm2calc_version_minor_match "${_gm2calc_version_header}")
  set(GM2Calc_VERSION_MINOR "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+GM2CALC_VERSION_PATCH[ \t]+([0-9]+)" _gm2calc_version_patch_match "${_gm2calc_version_header}")
  set(GM2Calc_VERSION_RELEASE "${CMAKE_MATCH_1}")

  set(GM2Calc_VERSION ${GM2Calc_VERSION_MAJOR}.${GM2Calc_VERSION_MINOR}.${GM2Calc_VERSION_RELEASE})

  if(GM2Calc_FIND_VERSION)
    if(GM2Calc_FIND_VERSION_EXACT AND NOT ${GM2Calc_VERSION} VERSION_EQUAL ${GM2Calc_FIND_VERSION})
      message(FATAL_ERROR "GM2Calc version ${GM2Calc_VERSION} found in ${GM2Calc_INCLUDE_DIRS}, "
        "but exact version ${GM2Calc_FIND_VERSION} is required.")
    elseif(${GM2Calc_VERSION} VERSION_LESS ${GM2Calc_FIND_VERSION})
      message(FATAL_ERROR "GM2Calc version ${GM2Calc_VERSION} found in ${GM2Calc_INCLUDE_DIRS}, "
        "but at least version ${GM2Calc_FIND_VERSION} is required.")
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GM2Calc
  FOUND_VAR GM2Calc_FOUND
  REQUIRED_VARS
    GM2Calc_LIBRARIES
    GM2Calc_INCLUDE_DIRS
)

if(GM2Calc_FOUND AND NOT TARGET GM2Calc::GM2Calc)
  add_library(GM2Calc::GM2Calc UNKNOWN IMPORTED)
  set_target_properties(GM2Calc::GM2Calc PROPERTIES
    IMPORTED_LOCATION "${GM2Calc_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${GM2Calc_INCLUDE_DIRS}"
  )
endif()

mark_as_advanced(
  GM2Calc_INCLUDE_DIRS
  GM2Calc_LIBRARIES
)
