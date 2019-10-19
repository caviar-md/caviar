
## WHAT IS THIS?

A post processing code for the output of 'plt_dealii' force_field.

This is a simple code to help calculating differential capacitance according to
equation (2) of the paper
dx.doi.org/10.1021/jp503224w | J. Phys. Chem. C 2014, 118, 18291−18298
by Merlet et. al.

which is,
  C = partial<q> / Kb T = ( <q^2> - <q>^2 ) / Kb T

a simple calculation!
  partial<q> = < (q - <q>) ^ 2) = < q^2 + <q>^2 - 2 q <q> > 
             = <q^2> + <q>^2 -2 <q>^2 = <q^2> - <q>^2

The files starting with 'o_' are some example data, input and output.

To compile,

 $ ./compile.sh

Then run,

 $ ./dif_capac

Follow the comments in the 'main_diff_capac_calculator.cpp 'and you will know 
what to do. If you have any questions, ask 'm.biagooi .at. gmail.com'

