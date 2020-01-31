
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

PROJECT(CAVIAR)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)


set (CAVIAR_MAJOR_VERSION 1)
set (CAVIAR_MINOR_VERSION 0)
set (CAVIAR_PATCH_VERSION 1)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# =======================================
# ======= add modules to the path =======
# =======================================

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

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

print_stars()

include(cmake/print_final_configuration.cmake)
include(cmake/print_final_message_before_build.cmake)

# =========================================
# ======= configure cpp source file =======
# =========================================

configure_file (
  "${CMAKE_SOURCE_DIR}/include/caviar/utility/caviar_config.h.in"
  "${CMAKE_SOURCE_DIR}/include/caviar/utility/caviar_config.h"
  )

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
endif (CAVIAR_WITH_DEALII)



# ========================================
# ========================================
# ========================================
