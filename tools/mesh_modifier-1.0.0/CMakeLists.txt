
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(MeshModifier)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)



FILE(GLOB _ALL_FILES include include/*)

FOREACH(_FILE ${_ALL_FILES})
  IF(IS_DIRECTORY ${_FILE})
    INCLUDE_DIRECTORIES(${_FILE})
  ENDIF()
ENDFOREACH()


FILE(GLOB_RECURSE sources 
     src/*.cpp 
     src/mesh_modifier/*.cpp)


ADD_LIBRARY(my_libs
  ${sources}
  )

  
ADD_EXECUTABLE(MeshModifier ./main.cpp)


TARGET_LINK_LIBRARIES(MeshModifier my_libs)
