
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

#include "caviar/objects/unique/molecule.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/unique/atom.h"
#include "caviar/interpreter/object_handler/preprocessors_new.h"
#include "caviar/objects/unique/molecule_group.h"


namespace caviar {
namespace objects {
namespace unique {


Molecule::Molecule (CAVIAR *fptr) : Unique{fptr},
    part_of_a_molecule_group{false}, upper_level_molecule_group{nullptr},
    position{Vector<double>{0,0,0}},
    velocity{Vector<double>{0,0,0}}
{
  // (maybe) the only user accessable creation funciton.
  FC_OBJECT_INITIALIZE_INFO
}    
       
   
Molecule::Molecule (const Molecule & m) : Unique{m},
atomic_bond{m.atomic_bond},
atomic_bond_index{m.atomic_bond_index},
atomic_angle{m.atomic_angle},
atomic_angle_index{m.atomic_angle_index}
{
  position = m.position;
  velocity = m.velocity;
  part_of_a_molecule_group = m.part_of_a_molecule_group;
  upper_level_molecule_group = m.upper_level_molecule_group;

  for (auto&& a : m.atoms) {
    Atom a_new (a);
    a_new.part_of_a_molecule = true; 

    // if not do this, 'a_new' parent molecule will be 'm' and not 'this'.
    a_new.upper_level_molecule = this;      
    atoms.push_back (a_new);
  }
}

Molecule::~Molecule () {}

void Molecule::verify_settings () {
  
}

bool Molecule::read ( caviar::interpreter::Parser * parser) {
  FC_OBJECT_READ_INFO

  while(true) {
    FC_IF_RAW_TOKEN_EOF_EOL
    FC_IF_GET_REAL3D(position)
    else if (string_cmp(ts,"add_atomic_bond") || string_cmp(ts,"atomic_bond") ) {
      objects::atom_data::Bond b;
      GET_OR_CHOOSE_A_INT(b.index_1,"","")
      GET_OR_CHOOSE_A_INT(b.index_2,"","")
      GET_OR_CHOOSE_A_INT(b.type,"","")
      GET_OR_CHOOSE_A_REAL(b.length,"","")

      atomic_bond.push_back(b);


      if (b.index_1 == b.index_2)
        error->all (FC_FILE_LINE_FUNC_PARSE, "bond indices cannot be similar.");

      bool index_1_exist = false;
      bool index_2_exist = false;
      for (unsigned int i = 0; i < atomic_bond_index.size(); ++i) {
        if (atomic_bond_index[i]==b.index_1) index_1_exist = true;
        if (atomic_bond_index[i]==b.index_2) index_2_exist = true;
      }
      if (!index_1_exist) atomic_bond_index.push_back(b.index_1);
      if (!index_2_exist) atomic_bond_index.push_back(b.index_2);


    } else if (string_cmp(ts,"add_atomic_angle") || string_cmp(ts,"atomic_angle") ) {
      objects::atom_data::Angle b;
      GET_OR_CHOOSE_A_INT(b.index_1,"","")
      GET_OR_CHOOSE_A_INT(b.index_2,"","")
      GET_OR_CHOOSE_A_INT(b.index_3,"","")
      GET_OR_CHOOSE_A_INT(b.type,"","")
      GET_OR_CHOOSE_A_REAL(b.value,"","")

      atomic_angle.push_back(b);


      if (b.index_1 == b.index_2 || b.index_1 == b.index_3 || b.index_2 == b.index_3)
        error->all (FC_FILE_LINE_FUNC_PARSE, "angle indices cannot be similar.");

      bool index_1_exist = false;
      bool index_2_exist = false;
      bool index_3_exist = false;
      for (unsigned int i = 0; i < atomic_angle_index.size(); ++i) {
        if (atomic_angle_index[i]==b.index_1) index_1_exist = true;
        if (atomic_angle_index[i]==b.index_2) index_2_exist = true;
        if (atomic_angle_index[i]==b.index_3) index_3_exist = true;
      }
      if (!index_1_exist) atomic_angle_index.push_back(b.index_1);
      if (!index_2_exist) atomic_angle_index.push_back(b.index_2);
      if (!index_3_exist) atomic_angle_index.push_back(b.index_3);


    } else if (string_cmp(ts,"output_xyz")) {
      output_xyz();
      continue;
    } else if (string_cmp(ts,"add_atom")) {
      add_atom (parser);
      return true;
    } else FC_ERR_UNDEFINED_VAR(ts)
  }
  
  return true;
}

bool Molecule::add_atom ( caviar::interpreter::Parser * parser) {

  FIND_OBJECT_BY_NAME(unique,it)
  FC_CHECK_OBJECT_CLASS_NAME(unique,it,atom)
  auto a =  *dynamic_cast<objects::unique::Atom *>(object_container->unique[it->second.index]);

  bool in_file = true;
  Vector<double> pos_ = {0,0,0};
  Vector<double> vel_ = {0,0,0};
  bool position_called = false;
  bool velocity_called = false;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"at_position")) {
      position_called = true;
      GET_OR_CHOOSE_A_REAL(pos_.x,"","")
      GET_OR_CHOOSE_A_REAL(pos_.y,"","")
      GET_OR_CHOOSE_A_REAL(pos_.z,"","")
    } else if (string_cmp(t,"at_velocity")) {
      velocity_called = true;
      GET_OR_CHOOSE_A_REAL(vel_.x,"","")
      GET_OR_CHOOSE_A_REAL(vel_.y,"","")
      GET_OR_CHOOSE_A_REAL(vel_.z,"","")
    } else FC_ERR_UNDEFINED_VAR(t)
  }


  if (position_called) a.position = pos_;
  if (velocity_called) a.velocity = vel_;

  a.part_of_a_molecule = true;
  a.upper_level_molecule = this;
  atoms.push_back (a);

  return in_file;
}


/*
void Molecule::correct_heritage () {
  for (unsigned int i = 0; i < atoms.size(); ++i) {
    atoms[i].upper_level_molecule = this;
    atoms[i].part_of_a_molecule = true; 
  }
}
*/

void Molecule::extract_all_e_pos_vel (std::vector<int>& e, std::vector<Vector<double>>&p,
 std::vector<Vector<double>>&v) {

  for (unsigned int i = 0; i < atoms.size(); ++i) {
    atoms[i].extract_all_e_pos_vel(e, p, v);
  }
}
 

void Molecule::output_xyz () {
  std::string out_file_name = "unique_molecule_" + object_name + ".xyz";
  std::ofstream out_file(out_file_name.c_str());
  out_file << atoms.size() << "\n" << "Atom\n";
  for (unsigned int i = 0; i < atoms.size(); ++i) {
    atoms[i].output_xyz (out_file);
  }
}

void Molecule::output_xyz (const std::string &st) {
  std::ofstream out_file(st.c_str());
  out_file << atoms.size() << "\n" << "Atom\n";
  for (unsigned int i = 0; i < atoms.size(); ++i) {
    atoms[i].output_xyz (out_file);
  }
}

void Molecule::output_xyz (std::ofstream & out_file) {
  for (unsigned int i = 0; i < atoms.size(); ++i) {
    atoms[i].output_xyz (out_file);
  }
}





bool Molecule::add_atom (const objects::unique::Atom &a, const Vector<double> &p, const Vector<double> &v) {
  unsigned int i = atoms.size();
  atoms.push_back (a);
  atoms[i].upper_level_molecule = this;
  atoms[i].part_of_a_molecule = true;
  atoms[i].position= p;
  atoms[i].velocity = v;
  return true; //WARNING
}


bool Molecule::add_atom (const objects::unique::Atom & a){
  unsigned int i = atoms.size();
  atoms.push_back (a);
  atoms[i].upper_level_molecule = this;
  atoms[i].part_of_a_molecule = true;
  return true; //WARNING  
}



Vector<double> Molecule::pos_tot () const {
  if (part_of_a_molecule_group) return position + upper_level_molecule_group->pos_tot();
   else return position;   
}

Vector<double> Molecule::vel_tot () const {
  if (part_of_a_molecule_group) return velocity + upper_level_molecule_group->vel_tot();
  else return velocity;  
}

} //unique
} //objects

} // namespace caviar

