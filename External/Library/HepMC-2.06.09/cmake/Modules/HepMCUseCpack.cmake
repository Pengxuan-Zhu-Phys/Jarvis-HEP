# define the environment for cpack
#


#
# Runs compiler with "-dumpversion" and parses major/minor
# version with a regex.
#
FUNCTION(_My_COMPILER_DUMPVERSION _OUTPUT_VERSION)

  EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
    ARGS ${CMAKE_CXX_COMPILER_ARG1} -dumpversion
    OUTPUT_VARIABLE _my_COMPILER_VERSION
  )
  set( COMPILER_VERSION ${_my_COMPILER_VERSION} PARENT_SCOPE)
  STRING(REGEX REPLACE "([0-9])\\.([0-9])(\\.[0-9])?" "\\1\\2"
    _my_COMPILER_VERSION ${_my_COMPILER_VERSION})

  SET(${_OUTPUT_VERSION} ${_my_COMPILER_VERSION} PARENT_SCOPE)
ENDFUNCTION()

#
# End functions/macros
#
#-------------------------------------------------------------------------------


macro( hepmc_find_compiler )
  if (My_COMPILER)
      SET (CPack_COMPILER_STRING ${My_COMPILER})
      message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
                     "using user-specified My_COMPILER = ${CPack_COMPILER_STRING}")
  else(My_COMPILER)
    # Attempt to guess the compiler suffix
    # NOTE: this is not perfect yet, if you experience any issues
    # please report them and use the My_COMPILER variable
    # to work around the problems.
    if (MSVC90)
      SET (CPack_COMPILER_STRING "-vc90")
    elseif (MSVC80)
      SET (CPack_COMPILER_STRING "-vc80")
    elseif (MSVC71)
      SET (CPack_COMPILER_STRING "-vc71")
    elseif (MSVC70) # Good luck!
      SET (CPack_COMPILER_STRING "-vc7") # yes, this is correct
    elseif (MSVC60) # Good luck!
      SET (CPack_COMPILER_STRING "-vc6") # yes, this is correct
    elseif (BORLAND)
      SET (CPack_COMPILER_STRING "-bcb")
    elseif("${CMAKE_CXX_COMPILER}" MATCHES "icl"
        OR "${CMAKE_CXX_COMPILER}" MATCHES "icpc")
      if(WIN32)
        set (CPack_COMPILER_STRING "-iw")
      else()
        set (CPack_COMPILER_STRING "-il")
      endif()
    elseif (MINGW)
        _My_COMPILER_DUMPVERSION(CPack_COMPILER_STRING_VERSION)
        SET (CPack_COMPILER_STRING "-mgw${CPack_COMPILER_STRING_VERSION}")
    elseif (UNIX)
      if (CMAKE_COMPILER_IS_GNUCXX)
          _My_COMPILER_DUMPVERSION(CPack_COMPILER_STRING_VERSION)
          # Determine which version of GCC we have.
	  if(APPLE)
              SET (CPack_COMPILER_STRING "-xgcc${CPack_COMPILER_STRING_VERSION}")
	  else()
              SET (CPack_COMPILER_STRING "-gcc${CPack_COMPILER_STRING_VERSION}")
	  endif()
      endif (CMAKE_COMPILER_IS_GNUCXX)
    endif()
    message(STATUS "Using compiler ${CPack_COMPILER_STRING}")
  endif(My_COMPILER)
endmacro( hepmc_find_compiler )


# note that hepmc_parse_version is used to define VERSION_MAJOR, etc.
set( CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR} )
set( CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR} )
set( CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH} )

set( CPACK_INCLUDE_TOPLEVEL_DIRECTORY 0 )
set( CPACK_GENERATOR TGZ )

hepmc_find_compiler()

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
  set(SLTYPE "slc")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set(SLTYPE "mac106")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
  set(SLTYPE "winxp")
endif ()
message(STATUS "CMAKE_HOST_SYSTEM_VERSION: ${CMAKE_HOST_SYSTEM_VERSION}")
message(STATUS "CMAKE_SYSTEM_VERSION: ${CMAKE_SYSTEM_VERSION}")
  if ( NOT CPack_COMPILER_STRING )
    set( PACKAGE_BASENAME ${CMAKE_SYSTEM_PROCESSOR}-${SLTYPE} )
  else ()
    set( PACKAGE_BASENAME ${CMAKE_SYSTEM_PROCESSOR}-${SLTYPE}${CPack_COMPILER_STRING} )
  endif ()
if ( NOT qualifier )
  set( CPACK_SYSTEM_NAME ${PACKAGE_BASENAME} )
else ()
  set( CPACK_SYSTEM_NAME ${PACKAGE_BASENAME}-${qualifier} )
endif ()
# check for extra qualifiers
if( NOT  CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE_TOLOWER default )
else()
   STRING(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_TOLOWER)
   if( ${CMAKE_BUILD_TYPE_TOLOWER} MATCHES "debug")
      set(CPACK_SYSTEM_NAME ${CPACK_SYSTEM_NAME}-debug )
   elseif( ${CMAKE_BUILD_TYPE_TOLOWER} MATCHES "release")
      set(CPACK_SYSTEM_NAME ${CPACK_SYSTEM_NAME}-opt )
   elseif( ${CMAKE_BUILD_TYPE_TOLOWER} MATCHES "minsizerel")
      set(CPACK_SYSTEM_NAME ${CPACK_SYSTEM_NAME}-prof )
   endif()   
endif()

message(STATUS "CPACK_SYSTEM_NAME = ${CPACK_SYSTEM_NAME}" )

include(CPack)
