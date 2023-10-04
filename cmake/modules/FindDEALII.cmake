
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


# ===============================
# =======  FINDING DEALII =======
# ===============================

FIND_PACKAGE(deal.II 8.5.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR={/PATH/TO/deal.II} to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ELSE()
  add_compile_definitions(DEALII_VERSION_MAJOR=${deal.II_VERSION_MAJOR})
  add_compile_definitions(DEALII_VERSION_MINOR=${deal.II_VERSION_MINOR})
  MESSAGE(STATUS "deal.II Version = ${deal.II_VERSION_MAJOR}.${deal.II_VERSION_MINOR}")
ENDIF()
 
DEAL_II_INITIALIZE_CACHED_VARIABLES()

MESSAGE(STATUS "deal.II deal.II_FOUND =" ${deal.II_FOUND}) # Writes 1 or 0 (I guess)
MESSAGE(STATUS "DEAL_II_DIR =" ${DEAL_II_DIR})  # writes nothig
MESSAGE(STATUS "Environment(DEAL_II_DIR) =" $ENV{DEAL_II_DIR})  # writes nothig
MESSAGE(STATUS "deal.II_DIR =" ${deal.II_DIR})  # writes something such as /usr/share/cmake/deal.II
#MESSAGE(STATUS "deal.II_LIBRARIES =" ${deal.II_LIBRARIES})  # writes nothig
#MESSAGE(STATUS "deal.II_INCLUDE_DIRS =" ${deal.II_INCLUDE_DIRS})  # writes nothig
#MESSAGE(STATUS "deal.II_DIR =" ${deal.II_LIBRARIES})  # writes nothig


#MESSAGE(STATUS "")
#MESSAGE(STATUS "deal.II deal.II_FOUND =" ${deal.II_FOUND}) # Writes 1 or 0 (I guess)


#MESSAGE(STATUS "")
#MESSAGE(STATUS "DEAL_II_INCLUDE_DIRS:" ${DEAL_II_INCLUDE_DIRS}) # Writes so many directories without space

#MESSAGE(STATUS "")
#MESSAGE(STATUS "DEAL_II_LIBRARIES:" ${DEAL_II_LIBRARIES}) # Writes so many libraries without space
# ===============================
# =======                 =======
# ===============================
