# - Find the MULTEM library
#
# Usage:
#   find_package(MULTEM [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   MULTEM_FOUND               ... true if multem is found on the system
#   MULTEM_LIB                 ... full path to multem library
#   MULTEM_INCLUDES            ... multem include directory
#
# The following variables will be checked by the function
#   MULTEM_ROOT               ... if set, the libraries are exclusively searched
#                               under this path
#   MULTEM_INCLUDE_DIR        ... multem include directory
#

#If environment variable MULTEMDIR is specified, it has same effect as MULTEM_ROOT
if( NOT MULTEM_ROOT AND ENV{MULTEMDIR} )
  set( MULTEM_ROOT $ENV{MULTEMDIR} )
endif()

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT MULTEM_ROOT )
  pkg_check_modules( PKG_MULTEM QUIET "multem" )
endif()

##Check whether to search static or dynamic libs
#set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )
#set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX} )

if( MULTEM_ROOT )

  #find libs
  find_library(
    MULTEM_LIB
    NAMES "multem"
    PATHS ${MULTEM_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
  )

  #find includes
  find_path(
    MULTEM_INCLUDES
    NAMES "multem/multem.h"
    PATHS ${MULTEM_ROOT}
    PATH_SUFFIXES "include"
    NO_DEFAULT_PATH
  )

else()

  find_library(
    MULTEM_LIB
    NAMES "libmultem"
    PATHS ${PKG_MULTEM_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
  )

  find_path(
    MULTEM_INCLUDES
    NAMES "multem/multem.h"
    PATHS ${PKG_MULTEM_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
  )

endif( MULTEM_ROOT )

# set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MULTEM DEFAULT_MSG MULTEM_INCLUDES MULTEM_LIB)

mark_as_advanced(MULTEM_INCLUDES MULTEM_LIB)

