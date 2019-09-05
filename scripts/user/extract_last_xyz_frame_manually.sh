#!/bin/bash

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


# Let's consider 'o_xyz.xyz' contains many frames of a simulation.
# with this command, one can extract the last frame of it to use it in another
# simulation as starting position.
# NUMBER is number of atoms plus 2;
# before using this, make sure the last frame is completed.

$ tail -n NUMBER o_xyz.xyz >o_xyz_last.xyz

# then use this

