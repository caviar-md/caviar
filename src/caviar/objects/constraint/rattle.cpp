
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Implementation of Rattle algorithm is done by Ashkan Shahmoradi
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

#include "caviar/objects/constraint/rattle.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"
#include "caviar/utility/interpreter_io_headers.h"

namespace caviar {
namespace objects {
namespace constraint {

Rattle::Rattle (CAVIAR *fptr) : Constraint{fptr},
    domain{nullptr},
    dt{-1.0}, error_tolerance{1e-6},
    domain_dh{Vector<double>{0,0,0}},
    domain_bc{Vector<int>{0,0,0}},
    initialized{false}
 {
  FC_OBJECT_INITIALIZE_INFO
  constraint_type = Constraint_t::Rattle;
}

Rattle::~Rattle () {}

bool Rattle::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"dt") ) {
      GET_OR_CHOOSE_A_REAL(dt,"","")
    } else if (string_cmp(t,"error_tolerance") ) {
      GET_OR_CHOOSE_A_REAL(error_tolerance,"","")
    } else if (string_cmp(t,"set_domain") || string_cmp(t,"domain")) {
      FIND_OBJECT_BY_NAME(domain,it)
      domain = object_container->domain[it->second.index];
      domain_dh = 0.5*(domain->upper_global - domain->lower_global);
      domain_bc = domain->boundary_condition;
    } 
    else FC_ERR_UNDEFINED_VAR(t)


  }
  return in_file;
}

void Rattle::verify_settings () {
  FC_NULLPTR_CHECK(atom_data)
  atom_data->record_owned_position_old = true;
  FC_NULLPTR_CHECK(domain)
  if (dt <= 0.0) error->all (FC_FILE_LINE_FUNC, "dt have to be a positive number");      

}


void Rattle::apply_on_position (int64_t) {
  FC_OBJECT_VERIFY_SETTINGS


  auto &vel = atom_data -> owned.velocity;
  auto &pos = atom_data -> owned.position;
  auto &pos_old = atom_data -> owned.position_old;
  auto &atomic_bond_vector = atom_data -> owned.atomic_bond_vector;	

  for (unsigned int i=0; i<atomic_bond_vector.size(); i++) { 

    auto Nc = atomic_bond_vector[i].size();
    if (Nc==0) continue;
    std::vector<double> C(Nc,0);

    double sum_err = 1.0;

    while(sum_err>error_tolerance) {

          
      for (unsigned int j=0; j<atomic_bond_vector[i].size(); j++) { 

        auto k1 = atomic_bond_vector[i][j].index_1, k2 = atomic_bond_vector[i][j].index_2;

        auto d = atomic_bond_vector[i][j].length;

        auto mass_inv_k1 = atom_data->owned.mass_inv[atom_data->owned.type[k1]];
        auto mass_inv_k2 = atom_data->owned.mass_inv[atom_data->owned.type[k2]];


        auto dr = domain -> fix_distance (pos[k1] - pos[k2]);
                
        auto dr_old = domain -> fix_distance (pos_old[k1] - pos_old[k2]);
                
        auto lambda = (dr*dr - (d*d)) /  (2.0*(mass_inv_k1 + mass_inv_k2)*(dr * dr_old));



                
        pos[k1] -= mass_inv_k1 *  lambda * dr_old;

        dr_old = domain -> fix_distance (pos_old[k2] - pos_old[k1]); // Note (k2 - k1)

        pos[k2] -= mass_inv_k2 *  lambda * dr_old;
                
        dr = domain -> fix_distance (pos[k1] - pos[k2]);



        auto dv = vel[k1] - vel[k2];

        auto etha =  (dr * dv) / ((mass_inv_k1+mass_inv_k2)*d*d);



        vel[k1] -= mass_inv_k1 * etha * dr;

        vel[k2] += mass_inv_k2 * etha * dr; // Note the Plus sign

        }


      sum_err = 0.0;
      for (unsigned int j=0; j<atomic_bond_vector[i].size(); j++) { 
        int k1 = atomic_bond_vector[i][j].index_1, k2 = atomic_bond_vector[i][j].index_2;

        auto d = atomic_bond_vector[i][j].length;

        auto dr = domain -> fix_distance (pos[k1] - pos[k2]);

        auto r2 = dr * dr;

        C[j] = (r2 - d*d)/(2*d*d);

        sum_err += C[j];
      }
      sum_err = abs( sum_err );	

      }
    }
}


} //constraint
} //objects
} // namespace caviar

