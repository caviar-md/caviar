
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

# ===========================================
# ======= built type check and config =======
# ===========================================

MACRO(check_build_type)

  if (CAVIAR_DEBUG_VERSION)

    set (CAVIAR_DEBUG_VERSION ON)
    set (CMAKE_BUILD_TYPE Debug)
    add_definitions(-DCAVIAR_DEBUG_VERSION)

  else()

    set (CAVIAR_DEBUG_VERSION OFF)
    set (CMAKE_BUILD_TYPE Release)

  endif()

  message(STATUS "CAVIAR_DEBUG_VERSION : ${CAVIAR_DEBUG_VERSION} ")
  message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE} ")

ENDMACRO()

# ===============================
# =======                 =======
# ===============================
