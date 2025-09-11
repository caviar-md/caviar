
#=======================================================================
#
# Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
#
# This file is part of the CAVIAR package.
#
# The CAVIAR package is free software; you can use it, redistribute
# it, and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either
# version 3.0 of the License, or (at your option) any later version.
# The full text of the license can be found in the file LICENSE at
# the top level of the CAVIAR distribution.
#
#=======================================================================

# ========================================
# ======= top level CMakeLists.txt =======
# ========================================

#CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

PROJECT(CAVIAR)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)


#set (CAVIAR_MAJOR_VERSION 1)
#set (CAVIAR_MINOR_VERSION 0)
#set (CAVIAR_PATCH_VERSION 1)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-cast-function-type") # -Wcast-function-type comes in openmpi compilation
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# =======================================
# ======= add modules to the path =======
# =======================================

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# ======================================= 
# ===== set default catched variables === 
# =======================================

# Default to Release if no build type is set
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type specified. Defaulting to Release.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type" FORCE)
endif()

# Let users pick valid build types in ccmake/cmake-gui
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")

# Enable debug build of CAVIAR (adds extra checks, logging, assertions)
set(CAVIAR_DEBUG_VERSION OFF CACHE BOOL "Build CAVIAR with debugging enabled")

# Enable OpenMP parallelization in CAVIAR
set(CAVIAR_WITH_OPENMP OFF CACHE BOOL "Enable OpenMP parallelization")

# Enable MPI parallelization in CAVIAR
set(CAVIAR_WITH_MPI OFF CACHE BOOL "Enable MPI parallelization")

# Use a single MPI domain for molecular dynamics (experimental mode)
set(CAVIAR_SINGLE_MPI_MD_DOMAIN OFF CACHE BOOL "Use a single MPI domain for MD simulation")

# Enable support for deal.II finite element library
set(CAVIAR_WITH_DEALII OFF CACHE BOOL "Enable deal.II integration")

# Enable deal.II with MPI support
set(CAVIAR_WITH_DEALII_MPI OFF CACHE BOOL "Enable deal.II with MPI support")

# Enable Eigen linear algebra library support
set(CAVIAR_WITH_EIGEN OFF CACHE BOOL "Enable Eigen library support")

# Path to Eigen installation (if not found automatically)
set(EIGEN_DIR "" CACHE PATH "Path to Eigen installation")

# Enable muParser mathematical expression parser
set(CAVIAR_WITH_MUPARSER OFF CACHE BOOL "Enable muParser support")

# Path to muParser installation (if not found automatically)
set(MUPARSER_DIR "" CACHE PATH "Path to muParser installation root")

# Path to muParser include directory (if not found automatically)
set(MUPARSER_INCLUDE_DIR "" CACHE PATH "Path to muParser include directory")

# ===============================
# ======= Load all macros =======
# ===============================

message(" ")
message("Load all macros:")
FILE(GLOB _macro_files "${CMAKE_CURRENT_SOURCE_DIR}/cmake/macros/*.cmake")
FOREACH(_file ${_macro_files})
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

print_stars()

# ========================================
# ========================================
# ========================================

check_build_type()

print_stars()

# ===================================
# ======= Load all configures =======
# ===================================

include(cmake/configure/configure_mpi.cmake)
include(cmake/configure/configure_dealii.cmake)
include(cmake/configure/configure_eigen.cmake)
include(cmake/configure/configure_openmp.cmake)
include(cmake/configure/configure_muparser.cmake)

print_stars()

include(cmake/print_final_configuration.cmake)
include(cmake/print_final_message_before_build.cmake)

# =========================================
# ======= configure cpp source file =======
# =========================================

#configure_file (
#  "${CMAKE_SOURCE_DIR}/include/caviar/utility/caviar_config.h.in"
#  "${CMAKE_SOURCE_DIR}/include/caviar/utility/caviar_config.h"
#  )

# =================================
# ======= Find Source files =======
# =================================

ADD_SUBDIRECTORY(src)

FILE(GLOB _ALL_FILES include)
FOREACH(_FILE ${_ALL_FILES})
  IF(IS_DIRECTORY ${_FILE})
    #MESSAGE( STATUS ${_FILE}  )
    INCLUDE_DIRECTORIES(${_FILE})
  ENDIF()
ENDFOREACH()


# ======================================
# ======= executable and linking =======
# ======================================

ADD_EXECUTABLE(CAVIAR src/main_CAVIAR.cpp)


if (CAVIAR_WITH_DEALII)
  DEAL_II_SETUP_TARGET(CAVIAR)
  TARGET_LINK_LIBRARIES(CAVIAR CAVIAR_lib)
else ()
  TARGET_LINK_LIBRARIES(CAVIAR CAVIAR_lib)
  if (CAVIAR_WITH_MUPARSER)
    target_link_libraries(CAVIAR muparser)
  endif()
endif (CAVIAR_WITH_DEALII)


if (CAVIAR_WITH_OPENMP)
  target_link_libraries(CAVIAR OpenMP::OpenMP_CXX)
endif()

#if (CAVIAR_WITH_MUPARSER)
#  target_link_libraries(CAVIAR muparser)
#endif()

# ========================================
# ========================================
# ========================================
