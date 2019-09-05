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


# ======
# q = q* q^
# q^ = 1.602 * 1e-19 Col =  1.602 * 1e-13 uC (i.e. micro-coloumb)
# =======
# q/A coefficitent (in (uC/cm^2) units): (q* * c = q/A)
# (1.602 * 1e-13 uC * 1e^+16 cm^(-2)= 1602 (uC/cm^2 )
# all :
#============================================

BID1=130.956
BID2=96.0
BID3=180.853
BID4=96.0
BID5=97.044
BID6=144.0
BIDS=600.004 # sum 1:5
BIDALL=730.96 # sum 1:6

C_TIME=0.000006   # time conversion coef (timestep to (ps))
CCCC=64.08 # 1602/25

#============================================
gnuplot << EOF
set term pdfcairo dashed enhanced size 6.75in,5in font "Helvetica,18"
#set term ps size 1024, 768
#set terminal eps
#set output 'o_induced_charges_1.eps'
set output 'o_induced_charges_1.pdf'
fz(x)=0
#set term png
#set output '$PNG_OUTPUT_NAME'
set xlabel "time (ns)"
set ylabel "charge density ({/Symbol m}C/{cm}^2)"
set key right center
set key title "Boundary ID   " 
plot [:15][:]\
'o_induced_charge' u (\$1*$C_TIME):((\$3)*$CCCC/$BID1) w l lw 2 ti '1', \
'o_induced_charge' u (\$1*$C_TIME):((\$4)*$CCCC/$BID2) w l lw 2 ti '2', \
'o_induced_charge' u (\$1*$C_TIME):((\$5)*$CCCC/$BID3) w l lw 2 ti '3', \
'o_induced_charge' u (\$1*$C_TIME):((\$6)*$CCCC/$BID4) w l lw 2 ti '4', \
'o_induced_charge' u (\$1*$C_TIME):((\$7)*$CCCC/$BID5) w l lw 2 ti '5', \
'o_induced_charge' u (\$1*$C_TIME):((\$8)*$CCCC/$BID6) w l lw 2 ti '6'
EOF

 
