
set term postscript enhanced eps color dl 1 size 10in,5in font "Helvetica,22"


set xlabel "x (Angstrom)"
set ylabel "potential"

set key left top
#set key title "Potential type" 

set output 'o_potential.eps'

pot_initial = 100.0
pot_debay=pot_initial/exp(1.0)

p 'data.dat' u 1:2 w lp lw 2, pot_debay w l lw 1

