
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

# ========================================================
# =======  CHECK DEALII REALATED FLAG DEPENDENCIES =======
# ========================================================

if (CAVIAR_WITH_DEALII_MPI)
  if (NOT CAVIAR_WITH_DEALII)
    set(CAVIAR_WITH_DEALII ON)
    message (STATUS "'CAVIAR_WITH_DEALII' is set to 'ON' because of a requested "
                    "flag 'CAVIAR_WITH_DEALII_MPI'.")
    message ("")
  endif()

  if (NOT CAVIAR_WITH_MPI)
    set(CAVIAR_WITH_MPI ON)
    message (STATUS "'CAVIAR_WITH_MPI' is set to 'ON' because of a requested "
                    "flag 'CAVIAR_WITH_DEALII_MPI'. Now re-including 'configue_mpi.cmake'")
    message ("")
    include(cmake/configure/configure_mpi.cmake)
  endif()

endif()


# ===================================================
# =======  ADDIND DEALII REALATED DEFINITIONS =======
# ===================================================

if (CAVIAR_WITH_DEALII)
  set (CAVIAR_WITH_MUPARSER ON)

  set (CAVIAR_WITH_DEALII ON)
  find_package(DEALII)
  add_definitions (-DCAVIAR_WITH_DEALII)
else ()
  set (CAVIAR_WITH_DEALII OFF)
endif ()

print_flag_value(CAVIAR_WITH_DEALII)


if (CAVIAR_WITH_DEALII_MPI)
  set (CAVIAR_WITH_DEALII_MPI ON)
  add_definitions (-DCAVIAR_WITH_DEALII_MPI)
else ()
  set (CAVIAR_WITH_DEALII_MPI OFF)
endif ()

print_flag_value(CAVIAR_WITH_DEALII_MPI)

# ==================================================
# ======= CHECK DEALII CONFIGURATION FOR MPI =======
# ==================================================

if (CAVIAR_WITH_DEALII_MPI)
# Are all dependencies fulfilled for using trilinos or petsc
  IF(NOT (DEAL_II_WITH_PETSC OR DEAL_II_WITH_TRILINOS) OR NOT DEAL_II_WITH_P4EST OR DEAL_II_PETSC_WITH_COMPLEX) # keep in one line
    MESSAGE(STATUS "In DEALII_WITH_MPI=TRUE mode, there has to be other flags: ")
    MESSAGE(FATAL_ERROR "
Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
    DEAL_II_WITH_PETSC = ON
    DEAL_II_PETSC_WITH_COMPLEX = OFF
    DEAL_II_WITH_P4EST = ON
or
    DEAL_II_WITH_TRILINOS = ON
    DEAL_II_WITH_P4EST = ON
One or both of these combinations are OFF in your installation but at least one is required for this tutorial step."
    )
  ENDIF()
endif()

# ===============================
# =======                 =======
# ===============================
