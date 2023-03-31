
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

#include "caviar/objects/constraint/atom_molarity.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/unique/atom.h"
#include "caviar/objects/unique/molecule.h"
#include "caviar/utility/interpreter_io_headers.h"

#include <random>
#include <algorithm>

CAVIAR_NAMESPACE_OPEN

namespace constraint {

static inline int int_floor(double x) 
{ 
    return (int)(x+100000) - 100000; 
}

Atom_molarity::Atom_molarity (CAVIAR *fptr) : Constraint{fptr} {
  FC_OBJECT_INITIALIZE_INFO
  atom_type = 0; 
  minimum_limit = 0;
  calculation_box_low = Vector<double> {0.0,0.0,0.0};
  calculation_box_high = Vector<double> {0.0,0.0,0.0};
  creation_box_low = Vector<double> {0.0,0.0,0.0};
  creation_box_high = Vector<double> {0.0,0.0,0.0};
  creation_try = 1;
  creation_molecule = nullptr;
  creation_atom = nullptr;
  settings_verified = false;
  steps = 0; check_steps=1;
  minimum_set = false;
  maximum_set = false;
  constraint_type = Constraint_t::Atom_molarity;
}

Atom_molarity::~Atom_molarity () {}

bool Atom_molarity::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"atom_type")) {
      GET_OR_CHOOSE_A_INT(atom_type,"","")
      if (atom_type < 0) error->all (FC_FILE_LINE_FUNC_PARSE, "atom_type have to be non-negative."); 
    } else if (string_cmp(t,"check_steps")) {
      GET_OR_CHOOSE_A_INT(check_steps,"","")
      if (check_steps < 0) error->all (FC_FILE_LINE_FUNC_PARSE, "check_steps have to be non-negative."); 
    } else if (string_cmp(t,"creation_try")) {
      GET_OR_CHOOSE_A_INT(creation_try,"","")
      if (creation_try < 0) error->all (FC_FILE_LINE_FUNC_PARSE, "creation_try have to be non-negative."); 
    } else if (string_cmp(t,"minimum_limit")) {
      GET_OR_CHOOSE_A_INT(minimum_limit,"","")
      if (minimum_limit < 0) error->all (FC_FILE_LINE_FUNC_PARSE, "minimum_limit have to be non-negative."); 
      minimum_set = true;
    } else if (string_cmp(t,"maximum_limit")) {
      GET_OR_CHOOSE_A_INT(maximum_limit,"","")
      if (maximum_limit < 0) error->all (FC_FILE_LINE_FUNC_PARSE, "maximum_limit have to be non-negative."); 
      maximum_set = true;
    } else if (string_cmp(t,"calculation_box")) {
      GET_OR_CHOOSE_A_REAL(calculation_box_low.x,"","")
      GET_OR_CHOOSE_A_REAL(calculation_box_high.x,"","")
      GET_OR_CHOOSE_A_REAL(calculation_box_low.y,"","")
      GET_OR_CHOOSE_A_REAL(calculation_box_high.y,"","")
      GET_OR_CHOOSE_A_REAL(calculation_box_low.z,"","")
      GET_OR_CHOOSE_A_REAL(calculation_box_high.z,"","")
    } else if (string_cmp(t,"creation_box")) {
      GET_OR_CHOOSE_A_REAL(creation_box_low.x,"","")
      GET_OR_CHOOSE_A_REAL(creation_box_high.x,"","")
      GET_OR_CHOOSE_A_REAL(creation_box_low.y,"","")
      GET_OR_CHOOSE_A_REAL(creation_box_high.y,"","")
      GET_OR_CHOOSE_A_REAL(creation_box_low.z,"","")
      GET_OR_CHOOSE_A_REAL(creation_box_high.z,"","")
    } else if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"set_creation_molecule") || string_cmp(t,"creation_molecule")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,molecule)
      unique::Molecule *a = dynamic_cast<unique::Molecule *>(object_container->unique[it->second.index]);
      creation_molecule = a;
    } else if (string_cmp(t,"set_creation_atom") || string_cmp(t,"creation_atom")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,atom)
      unique::Atom *a = dynamic_cast<unique::Atom *>(object_container->unique[it->second.index]);
      creation_atom = a;
    } else {
      error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
    }
  }
  return in_file;
}

void Atom_molarity::verify_settings () {
  FC_NULLPTR_CHECK(atom_data)

  if (creation_box_high.x < creation_box_low.x) 
    error->all (FC_FILE_LINE_FUNC, "creation_box_high.x < creation_box_low.x");
  if (creation_box_high.y < creation_box_low.y) 
    error->all (FC_FILE_LINE_FUNC, "creation_box_high.x < creation_box_low.x");
  if (creation_box_high.z < creation_box_low.z) 
    error->all (FC_FILE_LINE_FUNC, "creation_box_high.x < creation_box_low.x");

  if (calculation_box_high.x < calculation_box_low.x) 
    error->all (FC_FILE_LINE_FUNC, "calculation_box_high.x < calculation_box_low.x");
  if (calculation_box_high.y < calculation_box_low.y) 
    error->all (FC_FILE_LINE_FUNC, "calculation_box_high.x < calculation_box_low.x");
  if (calculation_box_high.z < calculation_box_low.z) 
    error->all (FC_FILE_LINE_FUNC, "calculation_box_high.x < calculation_box_low.x");

  if (creation_molecule==nullptr && creation_atom == nullptr)

    error->all (FC_FILE_LINE_FUNC, "creation_molecule==nullptr && creation_atom == nullptr"); 
  settings_verified = true;
}


// XXX not sure if this is the best place to implement the function
void Atom_molarity::apply (int64_t steps) { // III

  if ((steps%check_steps)!=0) { return;}


  // if nothing is set. Don't do any calculation.
  if (!(minimum_set || maximum_set)) return;

  if (!settings_verified) verify_settings();

  auto &pos = atom_data->owned.position;
  auto &type = atom_data->owned.type;


  // we record a list of atom indices that are in the area. It is used when we
  // have a maximum_limit constraint.
  std::vector <int> index_inside_domain;

  // for efficiency reasons
  if (maximum_set)
    index_inside_domain.reserve(maximum_limit+20);  


  // getting the number of atoms in the calculation_box
  auto dcal = calculation_box_high - calculation_box_low;
  int sum_of_type = 0;
  for (unsigned int i = 0; i < pos.size(); ++i) {
    auto pm = pos[i] - calculation_box_low;

    if (type[i] == static_cast<unsigned> (atom_type)) {
      if ((int_floor(pm.x / dcal.x) == 0) && 
          (int_floor(pm.y / dcal.y) == 0) &&
          (int_floor(pm.z / dcal.z) == 0) ) {
        ++sum_of_type;
        if (maximum_set) 
          index_inside_domain.push_back(i);
      }
    }
  }

  static std::mt19937 mt(1);

  // adding atoms if neccesary.
  if (minimum_set && sum_of_type < minimum_limit) {

  static std::uniform_real_distribution<double> dist_x(creation_box_low.x, creation_box_high.x);
  static std::uniform_real_distribution<double> dist_y(creation_box_low.y, creation_box_high.y);
  static std::uniform_real_distribution<double> dist_z(creation_box_low.z, creation_box_high.z);


  int tries = 0;

  if (creation_molecule!=nullptr) {
    auto a = *creation_molecule;
    while (tries < creation_try && sum_of_type < minimum_limit) {
      Vector<double> p {dist_x(mt),dist_y(mt), dist_z(mt)};
      a.position = p;
      if (atom_data->empty_of_atoms(a)) { 
        if (atom_data->add_molecule(a)) {
          sum_of_type++;
        }
      }
      ++tries;
    }
  } else if (creation_atom != nullptr) {
    auto a = *creation_atom;
    while (tries < creation_try && sum_of_type < minimum_limit) {
      Vector<double> p {dist_x(mt),dist_y(mt), dist_z(mt)};
      a.position = p;
      if (atom_data->empty_of_atoms(a)) { 
        if (atom_data->add_atom(a)) {
          sum_of_type++;
        }
      }
      ++tries;
    }
  } 
  }

  // removing atom if necessary
  if (maximum_set && sum_of_type > maximum_limit) {

    int diff = sum_of_type - maximum_limit;

    std::vector<int> v_delete_list;
    v_delete_list.reserve(diff);

    // Method I:
    // delete randomly. It may works, but it may be unnecessarily time consuming
    // in situations of low molarities. 
    /*
    std::uniform_int_distribution<> dist_i(0, index_inside_domain.size()-1);
    for (int n=0; n<diff; ++n) {
      int i = dist_i(mt);
      // element is in the vector
      if (std::find(v_delete_list.begin(), v_delete_list.end(), i) == v_delete_list.end())
        v_delete_list.push_back(index_inside_domain[i]);
    }*/
    

    // Method II:
    // Delete according to the index.
    // since the particles are identical we can easily remove them in a list.
    // Currently it makes no artifact in the simulations that I know.
    // we delete from the last to the first,
    // because in this case we will have less sorting steps. Since the 
    // 'index_inside_domain' vector is constructed from 0 to last;
    auto i_last = index_inside_domain.size() - 1;
    for (int n=0; n<diff; ++n) {
      v_delete_list.push_back(index_inside_domain[i_last - n]);
    }

    atom_data->remove_atom(v_delete_list);
    
  }

}

} //constraint

CAVIAR_NAMESPACE_CLOSE

