
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

#include "caviar/objects/atom_data.h"
#include "caviar/interpreter/communicator.h"
#include "caviar/interpreter/error.h"
#include "caviar/objects/unique/atom.h"
#include "caviar/objects/unique/atom_group.h"
#include "caviar/objects/unique/atom_list.h"
#include "caviar/objects/unique/molecule.h"
#include "caviar/objects/unique/molecule_group.h"
#include "caviar/objects/unique/molecule_list.h"

#include <algorithm>
#include <random>

namespace caviar {
namespace objects {



bool Atom_data::add_atom(caviar::objects::unique::Atom &a) {

  const auto id = get_global_id();
  const auto t = a.type;
  const auto pos = a.pos_tot ();
  const auto vel = a.vel_tot ();
  return add_atom (id, t, pos, vel);    
}

bool Atom_data::add_atom(caviar::objects::unique::Atom_group &ag) {
  for (auto&& a : ag.atoms)
    add_atom(a);
  return true; 
}

bool Atom_data::add_atom(caviar::objects::unique::Atom_list &al) {
  for (auto&& a : al.atoms)
    add_atom(*a);
  return true; 
}




bool Atom_data::add_molecule(caviar::objects::unique::Molecule &m) {

  std::vector<int> types;
  std::vector<Vector<double>> pos, vel;


  m.extract_all_e_pos_vel (types, pos, vel);

  // check all of atoms position validity...
  for (unsigned int i = 0; i < pos.size(); ++i) {
    if (!position_inside_local_domain(pos[i])) {
      return false;
    }
  }

  // adding the atoms while recording the map of molecule's atom index to atom_data index
  std::vector <int> indices (pos.size(), -1);
  for (unsigned int i = 0; i < pos.size(); ++i) {
    const auto id = get_global_id();
    const auto index = owned.position.size();    
    add_atom (id, types[i], pos[i], vel[i]);
    indices[i] = index;
  }


  // add the bonds if there are some. Changes the local indices to global ones
  // before addition.
  // not sure about this part to the end. But if it works, let it be!
  auto dummy_atomic_bond = m.atomic_bond;
  for (unsigned int j = 0;j<dummy_atomic_bond.size(); ++j) {
    dummy_atomic_bond[j].index_1 = indices[dummy_atomic_bond[j].index_1];
    dummy_atomic_bond[j].index_2 = indices[dummy_atomic_bond[j].index_2];
  }

  auto dummy_atomic_bond_index = m.atomic_bond_index;    

  for (unsigned int j = 0;j<dummy_atomic_bond_index.size(); ++j) {
    dummy_atomic_bond_index[j] = indices[dummy_atomic_bond_index[j]];
  }

  owned.atomic_bond_vector.push_back(dummy_atomic_bond); // XXX add condition of no zero sized push_back
  owned.atomic_bond_index_vector.push_back(dummy_atomic_bond_index);


  // add the angles if there are some. Changes the local indices to global ones
  // before addition.
  // not sure about this part to the end. But if it works, let it be!
  auto dummy_atomic_angle = m.atomic_angle;
  for (unsigned int j = 0;j<dummy_atomic_angle.size(); ++j) {
    dummy_atomic_angle[j].index_1 = indices[dummy_atomic_angle[j].index_1];
    dummy_atomic_angle[j].index_2 = indices[dummy_atomic_angle[j].index_2];
    dummy_atomic_angle[j].index_3 = indices[dummy_atomic_angle[j].index_3];
  }

  auto dummy_atomic_angle_index = m.atomic_angle_index;    

  for (unsigned int j = 0;j<dummy_atomic_angle_index.size(); ++j) {
    dummy_atomic_angle_index[j] = indices[dummy_atomic_angle_index[j]];
  }

  owned.atomic_angle_vector.push_back(dummy_atomic_angle);
  owned.atomic_angle_index_vector.push_back(dummy_atomic_angle_index);

  return true;
}

bool Atom_data::add_molecule(caviar::objects::unique::Molecule_group &mg) {
  for (auto&& a : mg.molecules)
    add_molecule(a);
  return true;
}

bool Atom_data::add_molecule(caviar::objects::unique::Molecule_list &ml) {
  for (auto&& a : ml.molecules)
    add_molecule(*a);
  return true;
}



} //objects

} // namespace caviar


