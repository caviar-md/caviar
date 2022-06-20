
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

# ===================================================
# ======= print_final_caviar_configuration =======
# ===================================================


print_stars()
  
message("Final CAVIAR configuration:")
message("")
print_flag_value(CAVIAR_DEBUG_VERSION)
print_flag_value(CMAKE_BUILD_TYPE)
message("")
print_flag_value(CAVIAR_WITH_OPENMP)
message("")
print_flag_value(CAVIAR_WITH_MPI)
print_flag_value(CAVIAR_SINGLE_MPI_MD_DOMAIN)
message("")
print_flag_value(CAVIAR_WITH_DEALII)
print_flag_value(CAVIAR_WITH_DEALII_MPI)
message("")
print_flag_value(CAVIAR_WITH_EIGEN)
print_flag_value(EIGEN_DIR)
message("")
print_flag_value(CAVIAR_WITH_MUPARSER)
print_flag_value(MUPARSER_DIR)
print_flag_value(MUPARSER_LIBRARY)
print_flag_value(MUPARSER_INCLUDE_DIR)
message("")
# ===============================
# =======                 =======
# ===============================
