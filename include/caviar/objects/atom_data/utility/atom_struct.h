
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// This file is part of the CAVIAR package.
//
// The CAVIAR package is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the CAVIAR distribution.
//
//========================================================================

#ifndef CAVIAR_OBJECTS_ATOMDATA_UTILITY_ATOMSTRUCT_H
#define CAVIAR_OBJECTS_ATOMDATA_UTILITY_ATOMSTRUCT_H

#include "caviar/utility/objects_common_headers.h"
#include "caviar/objects/atom_data/utility/bond.h"
#include "caviar/objects/atom_data/utility/angle.h"
#include "caviar/objects/atom_data/utility/proper_dihedral.h"

CAVIAR_NAMESPACE_OPEN

namespace atom_data
{
 /**
   * It contains all the physical data of atoms and molecules
   * but since it is used only once here, its name won't be of any use. Also
   * a good name would be 'atom_data' which is used before!.
   */
  struct Atom_struct  
  {
    /**
     * 'id' is a global and unique number assigned to an atom.
     */
    std::vector<GlobalID_t> id;

    /**
     * A tag to be done on the atoms.
     */
    std::vector<AtomType_t> tag;

    /**
     * Atom type decides the charge, mass and any other property shared between
     * a defined type (for example, Elements).
     */
    std::vector<AtomType_t> type;

    /**
     * 'mass' of an atom defined by the type. The mass may be used in
     * center-of-mass calculations and other functions. Do not depercate it.
     */
    std::vector<Real_t> mass;

    /**
     * simply the inverse value of 'mass' of an atom defined by the type.
     * since mass inverse is used in acceleration calculations.
     *
     */
    std::vector<Real_t> mass_inv;

    /**
     * 'charge' of an atom defined by the type.
     */
    std::vector<Real_t> charge;

    /**
     * 'radius' of an atom defined by the type. The user and the developers are
     * free to use this variable (for now!).
     */
    std::vector<Real_t> radius;

    /**
     * Different by atom_id. Can be changed in simulation (it is needed to be sent-recv. by MPI)
     */
    std::vector<Real_t> charge_atom;

    /**
     * Different by atom_id. Can be changed in simulation (it is needed to be sent-recv. by MPI)
     */
    std::vector<Real_t> mass_atom;

    /**
     * Atom kinematic properties in the current time-step.
     */
    std::vector<Vector<Real_t>> position, velocity, acceleration;

    /**
     * This vectors are used in some integrator schemes and constraint methods.
     * They can be defined in their related objects, but they may be needed in
     * more than one objects at once (for example, constraint::M_shake and
     * integrator::Leap_frog). This makes it the reason to define it here.
     * This function may be needed to have MPI_send-recv. process in these case.
     * look up to it.
     */
    std::vector<Vector<Real_t>> position_old, velocity_old, acceleration_old;

    /**
     * this vector is meaningful when there's one domain. We can calculate MSD
     * using this. It collects number of periodic domain cross for each particle.
     */
    std::vector<Vector<int>> msd_domain_cross;

    /**
     * this vector contain a molecule index for all the atoms. if it's '-1' the
     * atom is not of any molecule. This matters in the MPI process. All of the
     * atoms of a molecule should be existed in one process.
     */
    std::vector<int> molecule_index; //

    /**
     * number of total molecules.
     */
    int num_molecules;

    /**
     * The first std::vector, is the molecule index. the inner data contain bonds.
     */
    std::vector<std::vector<atom_data::Bond>> atomic_bond_vector;

    /**
     * Number of atomic bonds each atom have. It is used to limit
     * the bond creations.
     */
    std::vector<int> atomic_bond_count;

    /**
     * The first index, meaning
     * the first std::vector, is the molecule index. the inner data contain angles.
     */
    std::vector<std::vector<atom_data::Angle>> atomic_angle_vector;

    std::vector<std::vector<atom_data::Proper_dihedral>> atomic_properdihedral_vector;

  };
}
CAVIAR_NAMESPACE_CLOSE
#endif
