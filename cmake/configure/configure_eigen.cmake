
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


# ==============================
# =======  CONFIGURE EIGEN =======
# ==============================

if (CAVIAR_WITH_EIGEN)
  set (CAVIAR_WITH_EIGEN ON)

  find_package(EIGEN)

  add_definitions (-DCAVIAR_WITH_EIGEN)

else ()
  set (CAVIAR_WITH_EIGEN OFF)
endif ()

print_flag_value(CAVIAR_WITH_EIGEN)

if (CAVIAR_WITH_EIGEN)
  print_flag_value(EIGEN_DIR)
endif ()

# ===============================
# =======                 =======
# ===============================
