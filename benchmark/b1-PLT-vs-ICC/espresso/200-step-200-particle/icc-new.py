
# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

# THIS FILE IS GENERATED by modifications of an ESPResSo package script file.


import sys
import espressomd
import timeit

#from datetime import datetime


from espressomd.electrostatics import P3M
from espressomd.electrostatic_extensions import ICC

start_1 = timeit.default_timer()

S = espressomd.System(box_l=[1.0, 1.0, 1.0])
S.seed = 1
S.periodicity = [1, 1, 1]
# Parameters
box_l = 20.0
#nicc = 10
nicc = int(sys.argv[1])
q_test = 1.0
q_dist = 5.0

# System
S.box_l = [box_l, box_l, box_l + 5.0]
S.cell_system.skin = 0.4
S.time_step = 0.01

# ICC particles
nicc_per_electrode = nicc * nicc
nicc_tot = 2 * nicc_per_electrode
iccArea = box_l * box_l / nicc_per_electrode

iccNormals = []
iccAreas = []
iccSigmas = []
iccEpsilons = []

l = box_l / nicc
for xi in range(nicc):
    for yi in range(nicc):
        S.part.add(pos=[l * xi, l * yi, 0], q=-0.0001, fix=[1, 1, 1])
        iccNormals.append([0, 0, 1])

for xi in range(nicc):
    for yi in range(nicc):
        S.part.add(pos=[l * xi, l * yi, box_l],
                    q=0.0001, fix=[1, 1, 1])
        iccNormals.append([0, 0, -1])

iccAreas.extend([iccArea] * nicc_tot)
iccSigmas.extend([0] * nicc_tot)
iccEpsilons.extend([10000000] * nicc_tot)

# ADDING atoms
b2 = box_l * 0.5
num_of_atoms=9
l2 = (box_l - 2 ) / num_of_atoms
two_particle_test=1
if (two_particle_test==0) :
    S.part.add(pos=[b2, b2, b2 - q_dist / 2], q=q_test)
    S.part.add(pos=[b2, b2, b2 + q_dist / 2], q=-q_test)

else :
    for xi in range(num_of_atoms+1):
        for yi in range(num_of_atoms+1):
            S.part.add(pos=[1+l2 * xi, 1+l2 * yi, 12.5], q=-q_test)
            S.part.add(pos=[1+l2 * xi, 1+l2 * yi, 7.5], q=+q_test)



# Actors
#p3m = P3M(prefactor=1, mesh=32, cao=7, accuracy=1e-7)
p3m = P3M(prefactor=1, mesh=32, cao=7, accuracy=1e-5, alpha=1.112583061, r_cut=4.9, tune=False)
icc = ICC(
    n_icc=nicc_tot,
    convergence=1e-6,
    relaxation=0.75,
    ext_field=[0, 0, 0],
    max_iterations=100,
    first_id=0,
    eps_out=1,
    normals=iccNormals,
    areas=iccAreas,
    sigmas=iccSigmas,
    epsilons=iccEpsilons)

S.actors.add(p3m)
S.actors.add(icc)


#  time start


#start2=datetime.now()

#Statements





# Run
S.integrator.run()


fxyz= open("out.xyz","w+")

fxyz.write("%d\nAtom\n" %len(S.part))
for i in range(len(S.part)):
    fxyz.write("%d %f %f %f\n"  % (S.part[i].type, S.part[i].pos[0], S.part[i].pos[1], S.part[i].pos[2]))
fxyz.close()

# Analyze

fileq= open("bench.dat","a")

QL = sum(S.part[:nicc_per_electrode].q)
QR = sum(S.part[nicc_per_electrode:nicc_tot].q)
charge = (abs(QL) + abs(QR))/2.0

#start_2 = timeit.default_timer()

# Run
S.integrator.run(199)

#timestop
start_3 = timeit.default_timer()

fileq.write( "%d %d %.12f %.12f %.12f\n" % (nicc , len(S.part)-(2*nicc*nicc), charge, iccArea, start_3 - start_1) )

fileq.close()
#Your statements here
#print ("TT" ,  datetime.now()-start2)

#print('Time: ', stop - start)  

testcharge_dipole = q_test * q_dist
induced_dipole = 0.5 * (abs(QL) + abs(QR)) * box_l

# Result
#self.assertAlmostEqual(1, induced_dipole / testcharge_dipole, places=4)

# Test applying changes
#enegry_pre_change = S.analysis.energy()['total']
#pressure_pre_change = S.analysis.pressure()['total']
#icc.set_params(sigmas=[2.0] * nicc_tot)
#icc.set_params(epsilons=[20.0] * nicc_tot)
#enegry_post_change = S.analysis.energy()['total']
#pressure_post_change = S.analysis.pressure()['total']
#self.assertNotAlmostEqual(enegry_pre_change, enegry_post_change)
#self.assertNotAlmostEqual(pressure_pre_change, pressure_post_change)


