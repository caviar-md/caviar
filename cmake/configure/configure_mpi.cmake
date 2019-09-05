
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

if (CAVIAR_SINGLE_MPI_MD_DOMAIN)
  if (NOT CAVIAR_WITH_DEALII_MPI)
    set(CAVIAR_WITH_MPI 1)
    message (STATUS "'CAVIAR_WITH_MPI' is set to 'ON' because of a requested "
                    "flag 'CAVIAR_SINGLE_MPI_MD_DOMAIN'.")
    message ("")
  endif()
endif()

# ==============================
# =======  CONFIGURE MPI =======
# ==============================

if (CAVIAR_WITH_MPI)
  set (CAVIAR_WITH_MPI ON)

  find_package(MPI)

  add_definitions (-DCAVIAR_WITH_MPI)

  set (CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  #set (CMAKE_C_COMPILER ${MPI_C_COMPILER}) 

  message (STATUS "Set CXX_COMPILER TO ${MPI_CXX_COMPILER}")  
else ()
  set (CAVIAR_WITH_MPI OFF)
endif ()

print_flag_value(CAVIAR_WITH_MPI)

if (CAVIAR_SINGLE_MPI_MD_DOMAIN)
  set (CAVIAR_SINGLE_MPI_MD_DOMAIN ON)
  add_definitions (-DCAVIAR_SINGLE_MPI_MD_DOMAIN)
else()
  set (CAVIAR_SINGLE_MPI_MD_DOMAIN OFF)
endif()

print_flag_value(CAVIAR_SINGLE_MPI_MD_DOMAIN)

# ===============================
# =======                 =======
# ===============================
