
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
#include "caviar/utility/python_utils_def.h"
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


bool Atom::read ( caviar::interpreter::Parser *) {
  /*
  FC_OBJECT_READ_INFO
  while(true) {
    FC_IF_RAW_TOKEN_EOF_EOL
    FC_IF_GET_REAL3D(position)
    else FC_IF_GET_REAL3D(velocity)
    else FC_IF_GET_INT(type)
    else FC_ERR_UNDEFINED_VAR(ts)    
  }
  */
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

/*
FC_PYDEF_SETGET_PTR(Lj,atom_data,Atom_data);
FC_PYDEF_SETGET_PTR(Lj,domain,Domain);
FC_PYDEF_SETGET_PTR(Lj,neighborlist,Neighborlist);

FC_PYDEF_SETGET_STDVEC2D(Lj,epsilon,Real_t);  
FC_PYDEF_SETGET_STDVEC2D(Lj,sigma,Real_t);
FC_PYDEF_SETGET_STDVEC(Lj,epsilon_atom,Real_t);  
FC_PYDEF_SETGET_STDVEC(Lj,sigma_atom,Real_t);
FC_PYDEF_SETGET_STDVEC2D(Lj,cutoff_list,Real_t);
*/

FC_PYDEF_SETGET_CAVVEC(Atom,position,double);
FC_PYDEF_SETGET_CAVVEC(Atom,velocity,double);



void export_py_Atom () {

  using namespace boost::python;

  implicitly_convertible<std::shared_ptr<unique::Atom>,          
                         std::shared_ptr<Atom> >(); 

  class_<unique::Atom,boost::noncopyable>("Atom",init<caviar::CAVIAR*>())
  //class_<unique::Atom,boost::noncopyable>("Atom",init<caviar::CAVIAR*>())
  //class_<unique::Atom,std::shared_ptr<unique::Atom>>("Atom",init<caviar::CAVIAR*>())
  //class_<unique::Atom,A_Wrapper>("Atom",init<caviar::CAVIAR*>())
    .def_readwrite("type",&unique::Atom::type)    
    .add_property("position", &unique::Atom::get_position, &unique::Atom::set_position)
    .add_property("velocity", &unique::Atom::get_velocity, &unique::Atom::set_velocity)
  ;

  register_ptr_to_python<std::shared_ptr<unique::Atom> > ();
}



} //unique
} //objects

} // namespace caviar

