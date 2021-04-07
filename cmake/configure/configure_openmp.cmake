
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
# =======  CONFIGURE OPENMP                     =======
# =====================================================

if (CAVIAR_WITH_OPENMP)  
  set(CAVIAR_WITH_OPENMP ON)
  add_definitions (-DCAVIAR_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  message ("")  
  if (OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      #set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()  
else()
  set(CAVIAR_WITH_OPENMP OFF)
endif()



print_flag_value(CAVIAR_WITH_OPENMP)

# ===============================
# =======                 =======
# ===============================
