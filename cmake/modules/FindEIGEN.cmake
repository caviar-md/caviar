
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
# =======  FINDING EIGEN =======
# ============================

if (NOT EIGEN_DIR)
  SET( EIGEN_DIR "$ENV{EIGEN_DIR}" ) # Environment variable
ENDIF()

IF( NOT EIGEN_DIR )
  MESSAGE( FATAL_ERROR " Couldn't find Eigen library.
  Please point the environment variable EIGEN_DIR to the include directory of 
  your Eigen3 installation
or
  add the cmake definition with '-DEIGEN_DIR=THE_PATH/TO/EIGEN/LIBRARY' .
 ")
ENDIF()

if (NOT EXISTS "${EIGEN_DIR}/Eigen/Core")
  MESSAGE( FATAL_ERROR " Couldn't find Eigen Core header at: 
    '${EIGEN_DIR}/Eigen/Core' .")
ENDIF()

INCLUDE_DIRECTORIES ( "${EIGEN_DIR}" )

# ===============================
# =======                 =======
# ===============================
