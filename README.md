

## INTRODUCTION
CAVIAR is a molecular dynamics (MD) package developed by Morad Biagooi and 
Ehsan Nedaaee Oskoee. CAVIAR main goal is the simulation of soft materials at nano scale. However,
it can do much more. Handling complex geometries for MD is its speciallity.

The core is written in C++ with C++11 standard. The configurations for compile
is done using CMake. Some bash shell script are also used as tools for pre and
post processing. 

## CAVIAR SCRIPTING LANGUAGE (CASL)
CAVIAR has a built-in scripting language. It is one command per line. 
It is easy to learn and understand. See documentations and 
examples.


## HOW TO START?
The easiest way is to develop your simulation script from a working example.
Read the documentations and try to run your scripts accordingly.


## OS and LIBRARIES COMPATIBILITY
We have developed CAVIAR with g++ compiler on ubuntu OS:
  (g++ (Ubuntu 5.4.0-6ubuntu1~16.04.11) 5.4.0 20160609)
In addition to that, a deal.II (Version 8.5.1 ) library with for finite element 
calculations is used. 
Any lower version of these compilers or libraries may not work with CAVIAR,
specially if it does not support C++11 standards completely.
Higher versions may-or-maynot work depending on deprecations of the commands
that we have used. If you find any information about this, we gladly accept and
add it to our documentations.


## LICENSE
GNU Lesser General Public License as published by the Free Software Foundation.
Either version 3.0 or (at your option) any later version.
Note that:
THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.
EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER 
PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER 
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE 
QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE
DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

## Contact Information.
Morad Biagooi (m.biagooi .at. gmail.com)
Ehsan Nedaaee Oskoee (nedaaee .at. iasbs.ac.ir - nedaaee .at. gmail.com)



## Bug report
If you have encountered any problem while building or using CAVIAR, please
report it to 'm.biagooi .at. gmail.com'. We will have fix it as soon as possible.
There are three kind of bugs one can encounter:

1- Compile time bug: There should not have to be bug or error while building any
 realease version of CAVIAR. If something like that happened, there might be a problem
in compiler or library versions. 

2- Run time bug: If it is reported by an 'ERROR' massage of CAVIAR's
 interperter, please double check the script. If there's a segmentation fault,
you can use gdb to check the core file and find the problem. After that, contact
us and send your script to fix the problem.
Another type of bug is when the system freezes during an MPI run. 

3- Result bug: If your simulation result was not as it should be in a standard 
test simulation, there may be a problem in the parameters settings in the 
scripts. Also, there's a small possibility that the implemented algorithm is not
complete enough to support the physics of your problem. 




