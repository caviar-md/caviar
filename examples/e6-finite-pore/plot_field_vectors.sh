
set term postscript enhanced eps color dl 1 size 6.75in,5in font "Helvetica,22"


set xlabel "x (Angstrom)"
set ylabel "y (Angstrom)"

set key left top
#set key title "Potential type" 



set format x "%.f"

set xtic 5
set ytic 5
#set logscale y

ymax=1000

#2d
set output 'field_vectors_linear_filled2.eps'
plot [-10:10] [-5:5] "o_field_vectors_linear_filled2" u 1:2:($4*3):($5*3) w vectors ti 'smooth field'

set output 'field_vectors_linear_filled_si.eps'
plot [-10:10] [-5:5] "o_field_vectors_linear_filled_si" u 1:2:($4*3):($5*3) w vectors ti 'smooth field'

set output 'field_vectors_linear_filled_to.eps'
plot [-10:10] [-5:5] "o_field_vectors_linear_filled_to" u 1:2:($4*3):($5*3) w vectors ti 'smooth field'

set output 'field_vectors_linear_filled.eps'
plot [-31:24] [-12:12] "o_field_vectors_linear_filled" u 1:2:4:5 w vectors ti 'smooth field'

set output 'field_vectors_linear_empty.eps'
plot [-31:24] [-12:12] "o_field_vectors_linear_empty" u 1:2:4:5 w vectors ti 'smooth field'

set output 'field_vectors_linear_empty2.eps'
plot [-10:10] [-5:5] "o_field_vectors_linear_empty2" u 1:2:($4*3):($5*3) w vectors ti 'smooth field'

set output 'field_vectors_convex_filled.eps'
plot [-31:24] [-12:12] "o_field_vectors_convex_filled" u 1:2:4:5 w vectors ti 'smooth field'

set output 'field_vectors_linear_empty.eps'
plot [-31:24] [-12:12] "o_field_vectors_convex_empty" u 1:2:4:5 w vectors ti 'smooth field'

#3d
#splot [-31:24] [-12:12] "o_field_vectors" u 1:2:(0):4:5:(0) w vectors

