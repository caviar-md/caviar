
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

import caviarmd

c = caviarmd.caviar("")

#===== Atom definition =====

a1 = caviarmd.unique.Atom(c)
a1.type = 0
a1.position=[-1,0,0]

a2 = caviarmd.unique.Atom(c)
a2.type = 0
a2.position=[1,0,0]

#===== Domain

dom = caviarmd.domain.Box(c)
dom.lower_global = [-50,-50,-50]
dom.upper_global = [50,50,50]
dom.boundary_condition = [0,0,0]
dom.generate()

#========== Atom_data

adata = caviarmd.atom_data.Basic(c)

adata.ghost_cutoff = 5
adata.cutoff_extra = 0.01
adata.domain = dom
adata.add_atom(a1)
adata.add_atom(a2)
adata.add_type_mass(0,1.0)
adata.add_type_charge(0,0.0)

#===== Neighborlist

#neigh_verlet = caviarmd.neighborlist.Verlet_list(c)
#neigh_verlet.set_atom_data(adata)
#neigh_verlet.cutoff = 15
#neigh_verlet.dt = 0.001

#===== Neighborlist

#neigh_cell = caviarmd.neighborlist.Cell_list(c)
#neigh_cell.set_atom_data(adata)
#neigh_cell.cutoff=15
#neigh_cell.set_domain(dom)
#neigh_cell.make_neighlist = True
#neigh_cell.cutoff_neighlist=10
#===== force_field

#f_lj= caviarmd.force_field.Lj (c)

#f_lj.cutoff = 10.0
#f_lj.epsilon = [1.0]
#f_lj.sigma = [1.0]
#f_lj.set_neighborlist(neigh_verlet)
#f_lj.set_neighborlist(neigh_cell)
#f_lj.set_atom_data(adata)


#==== Integrator =====

#integ2 = caviarmd.integrator.Velocity_verlet  (c)
#integ2.set_atom_data(adata)
#integ2.dt=0.001

#====== writer 
#w1 = caviarmd.writer.Atom_data  (c)
#w1.set_atom_data (adata)
#w1.xyz_step = 200

#=====  simulator

#sim = caviarmd.md_simulator.Basic  (c)
#sim.set_integrator(integ2)
#sim.set_atom_data(adata)
#sim.add_force_field(f_lj)
#sim.add_neighborlist(neigh_verlet)
#sim.add_writer(w1)
#sim.initial_step = 0
#sim.final_step = 20000
#sim.dt = 0.001 
#sim.run()



