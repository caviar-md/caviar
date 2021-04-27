
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

# =====================================================
# =======  CHECK MPI REALATED FLAG DEPENDENCIES =======
# =====================================================


# =======================================
# =======  CONFIGURE BOOST PYTHON =======
# =======================================

if (CAVIAR_WITH_BOOST_PYTHON)
  set (CAVIAR_WITH_BOOST_PYTHON ON)

  
  #find_library(Bpython libboost_python36.so PATHS /usr/local/lib/)  
  
  #INCLUDE_DIRECTORIES(/usr/include/python3.6)
  

  FIND_PACKAGE(PythonInterp 3 REQUIRED)

  if (PYTHONINTERP_FOUND)
    if (UNIX AND NOT APPLE)
      if (PYTHON_VERSION_MAJOR EQUAL 3)
          FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX})
          FIND_PACKAGE(PythonInterp 3)
          FIND_PACKAGE(PythonLibs 3 REQUIRED)
      else()
          FIND_PACKAGE(Boost COMPONENTS python)
          FIND_PACKAGE(PythonInterp)
          FIND_PACKAGE(PythonLibs REQUIRED)
      endif()
    else()	
      if (PYTHON_VERSION_MAJOR EQUAL 3)
          FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
          FIND_PACKAGE(PythonInterp 3)
          FIND_PACKAGE(PythonLibs 3 REQUIRED)
      else()
          FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
          FIND_PACKAGE(PythonInterp)
          FIND_PACKAGE(PythonLibs REQUIRED)
      endif()
    endif()
  else()
      message("Python not found")
  endif()

  message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
  message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
  message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
  message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
  
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
  #LINK_LIBRARIES(${Boost_LIBRARIES} ${PYTHON_LIBRARIES}) # Deprecated but so convenient!

  add_definitions (-DCAVIAR_WITH_BOOST_PYTHON)

    
else ()
  set (CAVIAR_WITH_BOOST_PYTHON OFF)
endif ()

print_flag_value(CAVIAR_WITH_BOOST_PYTHON)


# ===============================
# =======                 =======
# ===============================
