
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

  find_package(BOOSTPYTHON)
  
  find_package(PYTHON)

  add_definitions (-DCAVIAR_WITH_BOOST_PYTHON)

  #set (CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  #set (CMAKE_C_COMPILER ${MPI_C_COMPILER}) 

  #message (STATUS "Set CXX_COMPILER TO ${MPI_CXX_COMPILER}")  
else ()
  set (CAVIAR_WITH_BOOST_PYTHON OFF)
endif ()

print_flag_value(CAVIAR_WITH_BOOST_PYTHON)


# ===============================
# =======                 =======
# ===============================
