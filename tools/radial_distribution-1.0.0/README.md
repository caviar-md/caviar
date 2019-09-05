
## WHAT IS THIS?

A post processing code for the output of xyz file formats.

It calculates radial distribution function for a system of atoms and molecules.


To compile,

 $ ./compile.sh

Then run,

 $ ./radial_distribution

It has some command line arguments.
If the input xyz file has velocity, the user have to set the configuration with
a '-v' argument:

 $ ./radial_distribution -v

It has other arguments. But the code flow is simple. Changing the C++ file is
straight forward.


