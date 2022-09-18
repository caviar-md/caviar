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

#================================================================
# Heatmap plotting example for potential exported from PLT force
# using gnuplot
#================================================================

set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
set output 'heatmaps-plot-potential.png'
set view map
set dgrid3d 200,200,2
set pm3d at b
set cbrange [-1000:1000]
set palette defined (-1000 "blue", 0 "white", 1000 "red")
unset key
unset surface
splot "o_potential_to" using 1:2:4
 
#================================================================
# the rest of the code are some useful commands, commented, from
# a good example on the internet downloaded from 
# 'http://gnuplot.sourceforge.net/demo/heatmaps.7.gnu'
#================================================================

# set format cb "%4.1f" 
# set view map scale 1
# set samples 25, 25
# set isosamples 50, 50
# set contour base
# set xyplane relative 0
# set cbtics border in scale 0,0 mirror norotate  autojustify
# set title "4D data (3D Heat Map)\nZ is contoured. Independent value is color-mapped" 
# set title  off set character 0, 1, 0 font "" textcolor lt -1 norotate
# set urange [ 5.00000 : 35.0000 ] noreverse nowriteback
# set vrange [ 5.00000 : 35.0000 ] noreverse nowriteback
# set xlabel "x" 
# set xlabel  off set character 3, 0, 0 font "" textcolor lt -1 norotate
# set xrange [ * : * ] noreverse writeback
# set x2range [ * : * ] noreverse writeback
# set ylabel "y" 
# set ylabel  off set character -1, 0, 0 font "" textcolor lt -1 norotate
# set yrange [ * : * ] noreverse writeback
# set y2range [ * : * ] noreverse writeback
# set zlabel "z" 
# set zlabel  off set character 2, 0, 0 font "" textcolor lt -1 norotate
# set zrange [ * : * ] noreverse writeback
# set cbrange [ * : * ] noreverse writeback
# set rrange [ * : * ] noreverse writeback
# set pm3d implicit at s
# set colorbox user
# set colorbox vertical origin screen 0.9, 0.2 size screen 0.03, 0.6 front  noinvert noborder
# sinc(x,y) = sin(sqrt((x-20.)**2+(y-20.)**2))/sqrt((x-20.)**2+(y-20.)**2)
# Z(x,y) = 100. * (sinc(x,y) + 1.5)
#color(x,y) = 10. * (1.1 + sin((x-20.)/5.)*cos((y-20.)/10.))
#NO_ANIMATION = 1
## Last datafile plotted: "++"
# splot '++' using 1:2:(Z($1,$2)):(color($1,$2)) with pm3d title "4 data columns x/y/z/color"

