
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


#include "caviar/objects/atom_data/basic.h"
#include "caviar/utility/python_utils_def.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/unique/atom.h"
#include "caviar/objects/unique/atom_group.h"
#include "caviar/objects/unique/atom_list.h"
#include "caviar/objects/unique/molecule.h"
#include "caviar/objects/unique/molecule_group.h"
#include "caviar/objects/unique/molecule_list.h"
#include "caviar/objects/neighborlist/cell_list.h"
#include "caviar/interpreter/object_handler/preprocessors_new.h"
#include "caviar/objects/domain.h"

namespace caviar {
namespace objects {
namespace atom_data {

Basic::Basic (CAVIAR *fptr) : Atom_data{fptr} {
  FC_OBJECT_INITIALIZE_INFO
}

void Basic::allocate () {
  
}
/*
#define FIND_UNIQUE_OBJECT_BY_NAME_NIC(OBJECT_TYPE,ITERATOR_NAME) \
  std::string name_to_find_;\
  MAKE_ERROR_MASSAGE_EXPECTED(err,OBJECT_TYPE,"name.")\
  GET_A_STRING(name_to_find_,"", err)\
  CHECK_NAME_EXISTANCE(name_to_find_, ITERATOR_NAME, "","")\
  if (ITERATOR_NAME->second.type != caviar::interpreter::object_handler::gdst( #OBJECT_TYPE ))\
    error->all(FC_FILE_LINE_FUNC_PARSE,": undefined object. ");

#define FIND_UNIQUE_OBJECT_BY_NAME(OBJECT_TYPE,ITERATOR_NAME) \
  std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator ITERATOR_NAME;\
  FIND_UNIQUE_OBJECT_BY_NAME_NIC(OBJECT_TYPE,ITERATOR_NAME)
*/
bool Basic::read (caviar::interpreter::Parser *) {
    /*
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;

    if (string_cmp(t,"ghost_cutoff")) {
      GET_OR_CHOOSE_A_REAL(ghost_cutoff,"","")
      if (ghost_cutoff < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "ghost_cutoff have to non-negative."); 
    } else if (string_cmp(t,"neighborlist_cutoff")) {
      GET_OR_CHOOSE_A_REAL(neighborlist_cutoff,"","")
      if (neighborlist_cutoff < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "neighborlist_cutoff have to non-negative."); 
    } else if (string_cmp(t,"cutoff_extra")) {
      GET_OR_CHOOSE_A_REAL(cutoff_extra,"","")
      if (cutoff_extra < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "neighborlist_cutoff have to non-negative."); 
    } else if (string_cmp(t,"k_b")) {
      GET_OR_CHOOSE_A_REAL(k_b,"","")
      if (k_b < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "k_b have to non-negative."); 
    } else if (string_cmp(t,"add_atom")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,atom)
      objects::unique::Atom *a = dynamic_cast<objects::unique::Atom *>(object_container->unique[it->second.index]);
      add_atom(*a);
    } else if (string_cmp(t,"add_atom_group")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,atom_group)
      objects::unique::Atom_group *a = dynamic_cast<objects::unique::Atom_group *>(object_container->unique[it->second.index]);
      add_atom(*a);
    } else if (string_cmp(t,"add_atom_list")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,atom_list)
      objects::unique::Atom_list *a = dynamic_cast<objects::unique::Atom_list *>(object_container->unique[it->second.index]);
      add_atom(*a);
    } else if (string_cmp(t,"add_molecule")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,molecule)
      objects::unique::Molecule *a = dynamic_cast<objects::unique::Molecule *>(object_container->unique[it->second.index]);
      add_molecule(*a);
    } else if (string_cmp(t,"add_molecule_group")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,molecule_group)
      objects::unique::Molecule_group *a = dynamic_cast<objects::unique::Molecule_group *>(object_container->unique[it->second.index]);
      add_molecule(*a);
    } else if (string_cmp(t,"add_molecule_list")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,molecule_list)
      objects::unique::Molecule_list *a = dynamic_cast<objects::unique::Molecule_list *>(object_container->unique[it->second.index]);
      add_molecule(*a);
    } else if (string_cmp(t,"add_type_radius")) {
      //auto ind = parser->get_positive_int();
      //auto r = parser->get_real();
      int ind = 0;
      double r = 0;
      GET_OR_CHOOSE_A_INT(ind,"","")
      GET_OR_CHOOSE_A_REAL(r,"","")
      if (static_cast<int> (owned.radius.size()) < ind+1) {
        owned.radius.resize(ind+1);
      }
      owned.radius[ind] = r;
      if (r<0) output->warning("you have entered a negative value for atom radius.");
    } else if (string_cmp(t,"add_type_charge")) {
      //auto ind = parser->get_positive_int();
      //auto c = parser->get_real();
      int ind = 0;
      double c = 0;
      GET_OR_CHOOSE_A_INT(ind,"","")
      GET_OR_CHOOSE_A_REAL(c,"","")
      if (static_cast<int> (owned.charge.size()) < ind+1) {
        owned.charge.resize(ind+1);
      }
      owned.charge[ind] = c;
    } else if (string_cmp(t,"add_type_mass")) {

      //auto ind = parser->get_positive_int();
      //auto m = parser->get_real();
      int ind = 0;
      double m = 0;
      GET_OR_CHOOSE_A_INT(ind,"","")
      GET_OR_CHOOSE_A_REAL(m,"","")
      if (static_cast<int> (owned.mass.size()) < ind+1) {
        owned.mass.resize(ind+1);
        owned.mass_inv.resize(ind+1);
      }
      owned.mass[ind] = m;
      if (m==0.0) owned.mass_inv[ind] = 0.0;
      else owned.mass_inv[ind] = 1.0/m;
      if (m<0) output->warning("you have entered a negative value for atom mass.");
    } else if (string_cmp(t,"set_domain") || string_cmp(t,"domain")) {
      FIND_OBJECT_BY_NAME(domain,it)
      domain = object_container->domain[it->second.index];
    } else if (string_cmp(t,"set_cell_list") || string_cmp(t,"cell_list")) {
      FIND_OBJECT_BY_NAME(neighborlist,it)
      FC_CHECK_OBJECT_CLASS_NAME(neighborlist,it,cell_list)
      cell_list = dynamic_cast<objects::neighborlist::Cell_list *>(object_container->neighborlist[it->second.index]);
    } else if (string_cmp(t,"add_xyz_data_file")) {
      add_xyz_data_file(parser);
      return true;
    } else if (string_cmp(t,"set_owned_position")) {
      auto ind = parser->get_int();
      auto x = parser->get_real();
      auto y = parser->get_real();
      auto z = parser->get_real();
      owned.position[ind] = Vector<Real_t> {x, y, z};
    } else if (string_cmp(t,"set_owned_velocity")) {
      auto ind = parser->get_int();
      auto x = parser->get_real();
      auto y = parser->get_real();
      auto z = parser->get_real();
      owned.velocity[ind] = Vector<Real_t> {x, y, z};
    } else if (string_cmp(t,"set_owned_acceleration")) {
      auto ind = parser->get_int();
      auto x = parser->get_real();
      auto y = parser->get_real();
      auto z = parser->get_real();
      owned.acceleration[ind] = Vector<Real_t> {x, y, z};
    } else if (string_cmp(t,"add_random_velocity")) {
      add_random_velocity();
      return true;
    } else if (string_cmp(t,"n_r_df")) {
      GET_OR_CHOOSE_A_INT(n_r_df,"","")
    } else FC_ERR_UNDEFINED_VAR(t)
   
  }
  return in_file;
 */  
  return true;
}

void Basic::output_data (int ) {
  
}

FC_PYDEF_SETGET_PTR(Basic,domain,Domain);
/*
FC_PYDEF_SETGET_PTR(Lj,atom_data,Atom_data);

FC_PYDEF_SETGET_PTR(Lj,neighborlist,Neighborlist);

FC_PYDEF_SETGET_STDVEC2D(Lj,epsilon,Real_t);  
FC_PYDEF_SETGET_STDVEC2D(Lj,sigma,Real_t);
FC_PYDEF_SETGET_STDVEC(Lj,epsilon_atom,Real_t);  
FC_PYDEF_SETGET_STDVEC(Lj,sigma_atom,Real_t);
FC_PYDEF_SETGET_STDVEC2D(Lj,cutoff_list,Real_t);

adata.ghost_cutoff = 5
adata.cutoff_extra = 0.01
adata.set_domain(dom)
adata.add_atom(a1)
adata.add_atom(a2)
adata.add_type_mass(0,1.0)
adata.add_type_charge(0,0.0)

*/

//virtual bool add_atom(const std::shared_ptr<caviar::objects::unique::Atom> &a);
//virtual bool add_atom(caviar::objects::unique::Atom_group &a);
//virtual bool add_atom(caviar::objects::unique::Atom_list &a);

// thin wrapper
//bool add_atom1(const std::shared_ptr<caviar::objects::unique::Atom> &a) {
//  return Basic::add_atom(a);
//}
// manual wrapper
//bool (Basic::*add_atom_ptr1) (const std::shared_ptr<caviar::objects::unique::Atom> &a) = &Basic::add_atom;

void export_py_Basic () {


  

  using namespace boost::python;

  implicitly_convertible<std::shared_ptr<atom_data::Basic>,          
                         std::shared_ptr<Atom_data> >(); 

  class_<atom_data::Basic,boost::noncopyable>("Basic",init<caviar::CAVIAR*>())
    .def_readwrite("ghost_cutoff",&atom_data::Basic::ghost_cutoff)    
    .def_readwrite("cutoff_extra",&atom_data::Basic::cutoff_extra)
    .def("add_atom",&atom_data::Basic::add_atom_wrap1)
    .def("add_type_mass",&atom_data::Basic::add_masses)
    .def("add_type_charge",&atom_data::Basic::add_charges)
    .add_property("domain",&atom_data::Basic::get_domain,&atom_data::Basic::set_domain)          
  ;

}


} //atom_data
} //objects
} // namespace caviar

