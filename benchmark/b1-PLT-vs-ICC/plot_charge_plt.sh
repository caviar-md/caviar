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


set term postscript enhanced eps color dl 4 size 6.75in,5in font "Helvetica,30"

set xtics 1

set xlabel "Number of Global Refinements"
set ylabel "Induced Charge (dimensionless)"
set key left top
set pointsize 2
set output 'cavi-induced-charge.eps'
plot [:][-1:]\
'caviar/cavi-200-step-200-particle.dat' u 1:3  w lp lw 2   pt 4 lt rgb "black" ti 'case (i)', \
'caviar/cavi-200-step-100-particle-random.dat' u 1:3  w lp lw 2 pt 6  lt rgb "red"  ti 'case (ii)', \


set xlabel "Number of Global Refinements"
set ylabel "CPU-time (second)"
set key left top
set pointsize 2
set output 'cavi-cpu-time.eps'
plot [:][:]\
'caviar/cavi-200-step-200-particle.dat' u 1:($2/200)  w lp lw 2  pt 4  lt rgb "black" ti 'case (i)', \
'caviar/cavi-200-step-100-particle-random.dat' u 1:($2/200)  w lp lw 2 pt 6   lt rgb "red" ti 'case (ii)', \

