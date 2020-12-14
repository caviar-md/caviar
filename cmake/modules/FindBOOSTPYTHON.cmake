
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

# =====================================
# =======  FINDING BOOST PYTHON =======
# =====================================

# to be able to call the system FindMPI.cmake module, we temporarily disable
# ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)


# =========================================
# ======= find boost python =======
# =========================================

find_package( Boost COMPONENTS python3 REQUIRED )

if( Boost_FOUND )
        MESSAGE("********************************FOUND BOOST***********************")
endif( Boost_FOUND )           


# NOT TESTED
set (Boost_VERBOSE ON)

FIND_PACKAGE(boost_python REQUIRED
  HINTS ${BOOST_PYTHON_DIR} $ENV{BOOST_PYTHON_DIR}
)

IF(NOT ${boost_python_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "Could not locate a (sufficiently recent) version of boost python\n\n"
    "You may want to either pass a flag -DBOOST_PYTHON_DIR=/path/to/boost_python to cmake\n"
    "or set an environment variable \"BOOST_PYTHON_DIR\" that contains this path."
    )    
endif()

message(STATUS "Boost_LIBRARY_DIRS :${Boost_LIBRARY_DIRS}")
message(STATUS "Boost_INCLUDE_DIRS : ${Boost_INCLUDE_DIRS}")

message(STATUS "BOOST_LIBRARIES : ${BOOST_LIBRARIES}")
message(STATUS "boost_python_DIR: ${boost_python_DIR}")

message(STATUS "boost_python: ${boost_python}")

message(STATUS boost_python_FOUND : ${boost_python_FOUND})
message(STATUS boost_python_INCLUDE_DIRS : ${boost_python_INCLUDE_DIRS})
message(STATUS boost_python_INCLUDES : ${boost_python_INCLUDES})
message(STATUS boost_python_LIBRARIES : ${boost_python_LIBRARIES})
message(STATUS boost_python_LIBS : ${boost_python_LIBS})
message(STATUS boost_python_DEFINITIONS : ${boost_python_DEFINITIONS})

#find_library(Bpython libboost_python36.so PATH /usr/local/lib) #boostpython

# now, we re-add the local module path
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# ===============================
# =======                 =======
# ===============================
