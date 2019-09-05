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


#plot 'o_mesh_boundary_normals' with the command below on gnuplot:

splot 'o_mesh_boundary_normals' u 1:2:3:4:5:6 w vectors

#The arrows has to point toward inside the simulation geometry (box or anything)
