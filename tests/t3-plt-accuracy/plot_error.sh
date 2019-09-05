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

#============================================
gnuplot << EOF
set term pdfcairo dashed enhanced size 6.75in,5in font "Helvetica,18"
set output 'o_error_f.pdf'
fz(x)=0
set xlabel "r/R"
set ylabel "Force error"
plot [:][0:0.5] \
'refine2.txt' u 1:15 w l lw 2 ti '2', \
'refine3.txt' u 1:15 w l lw 2 ti '3', \
'refine4.txt' u 1:15 w l lw 2 ti '4'
EOF


#============================================
gnuplot << EOF
set term pdfcairo dashed enhanced size 6.75in,5in font "Helvetica,18"
set output 'o_error_p.pdf'
fz(x)=0
#set term png
#set output '$PNG_OUTPUT_NAME'
set xlabel "r/R"
set ylabel "Potential error"
plot [:][0:0.5] \
'refine2.txt' u 1:14 w l lw 2 ti '2', \
'refine3.txt' u 1:14 w l lw 2 ti '3', \
'refine4.txt' u 1:14 w l lw 2 ti '4'
EOF

 
