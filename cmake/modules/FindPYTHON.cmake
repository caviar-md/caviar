
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
# =======  FINDING PYTHON =======
# =====================================

# to be able to call the system FindMPI.cmake module, we temporarily disable
# ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)


# =========================================
# ======= find python =====================
# =========================================

find_package(Python3 COMPONENTS Interpreter Development)

execute_process(COMMAND ${Python3_EXECUTABLE}
                -c "import distutils.sysconfig as cg; print(cg.get_python_inc())"
                OUTPUT_VARIABLE PYTHON_INCLUDE_PATH
                OUTPUT_STRIP_TRAILING_WHITESPACE)

set(PYTHON_INCLUDE_PATH ${Python3_INCLUDE_DIRS} CACHE PATH "Python Include Directory")
mark_as_advanced(PYTHON_INCLUDE_PATH)

message(STATUS "PYTHON_INCLUDE_PATH = ${PYTHON_INCLUDE_PATH}")
include_directories(${PYTHON_INCLUDE_PATH})

if(NOT PYTHON_INSTDIR)
    execute_process(COMMAND ${Python3_EXECUTABLE}
                -c "import distutils.sysconfig as cg; print(cg.get_python_lib(1,0,prefix='${CMAKE_INSTALL_EXEC_PREFIX}'))"
                OUTPUT_VARIABLE PYTHON_INSTDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

execute_process(COMMAND ${Python3_EXECUTABLE}
                        -c "import sys; print(sys.version[:3])"
                        OUTPUT_VARIABLE PYTHON_VERSION
                        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND ${Python3_EXECUTABLE}
                        -c "import sys; print(sys.version[:3].replace('.', ''))"
                        OUTPUT_VARIABLE PYTHON_VERSION_NO_DOT
                        OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "PYTHON_VERSION = ${PYTHON_VERSION}")
message(STATUS "PYTHON_LIBDIR = ${Python3_LIBRARY_DIRS}")

set(PYTHON_LIBRARIES ${Python3_LIBRARIES})

if(NOT Python3_LIBRARIES)
  message(FATAL_ERROR "Python libraries not found!")
endif(NOT Python3_LIBRARIES)
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")

mark_as_advanced(PYTHON_INCLUDE_PATH PYTHON_LIBRARIES)

# now, we re-add the local module path
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# ===============================
# =======                 =======
# ===============================
