
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

#include "caviar/objects/unique/atom.h"
#include "caviar/objects/unique/atom_group.h"
#include "caviar/objects/unique/molecule.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/interpreter/object_handler/preprocessors_new.h"

namespace caviar {
namespace objects {
namespace unique {

Atom::Atom (CAVIAR *fptr) : Unique{fptr},
    part_of_a_molecule{false}, upper_level_molecule{nullptr},
    part_of_a_atom_group{false}, upper_level_atom_group{nullptr},
    position{Vector<double>{0,0,0}},
    velocity{Vector<double>{0,0,0}}, type{0}  {
  // this is the user accessable creation function. So we only
  // call the macro initializer here.
  FC_OBJECT_INITIALIZE_INFO
}   
   
  
Atom::~Atom () {}

void Atom::verify_settings () {
  
}


Atom::Atom (const Atom & a) : Unique{a} {
  part_of_a_molecule = a.part_of_a_molecule;
  upper_level_molecule = a.upper_level_molecule;
  part_of_a_atom_group = a.part_of_a_atom_group;
  upper_level_atom_group = a.upper_level_atom_group;
  position = a.position; 
  velocity = a.velocity;
  type = a.type;
}


bool Atom::read ( caviar::interpreter::Parser * parser) {
  FC_OBJECT_READ_INFO
  while(true) {
    FC_IF_RAW_TOKEN_EOF_EOL
    FC_IF_GET_REAL3D(position)
    else FC_IF_GET_REAL3D(velocity)
    else FC_IF_GET_INT(type)
    else FC_ERR_UNDEFINED_VAR(ts)    
  }
  return true;
}

void Atom::extract_all_e_pos_vel (std::vector<int>& e, std::vector<Vector<double>>&p,
 std::vector<Vector<double>>&v) {
  e.push_back (type);
  p.push_back (pos_tot());
  v.push_back (vel_tot());
  //std::cout << "1pos : " << pos_tot() << std::endl;
}

void Atom::output_xyz (std::ofstream & out_file) {
  const auto p = pos_tot();
  out_file << type << " " << p.x << " " << p.y << " " << p.z << std::endl;  

}



Vector<double> Atom::pos_tot () const {
  if (part_of_a_molecule) return position + upper_level_molecule->pos_tot();
  else if (part_of_a_atom_group) return position + upper_level_atom_group->pos_tot();
  else return position;   
}

Vector<double> Atom::vel_tot () const {
  if (part_of_a_molecule) return velocity + upper_level_molecule->vel_tot();
  else if (part_of_a_atom_group) return velocity + upper_level_atom_group->vel_tot();
  else return velocity;    
}



} //unique
} //objects

} // namespace caviar

