
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

# ============================
# =======  FINDING MPI =======
# ============================

# to be able to call the system FindMPI.cmake module, we temporarily disable
# ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)


# Call the system FindMPI.cmake module:
find_package (MPI)

if (NOT MPI_CXX_FOUND)
  set (CAVIAR_WITH_MPI MPI-NOTFOUND)
endif ()

# now, we re-add the local module path
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# ===============================
# =======                 =======
# ===============================
