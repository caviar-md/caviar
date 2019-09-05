
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Implementation of M-Shake algorithm is done by Ashkan Shahmoradi
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

#include "caviar/objects/constraint/m_shake.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"
#include "caviar/utility/interpreter_io_headers.h"
#include <string>

namespace caviar {
namespace objects {
namespace constraint {



M_shake::M_shake (CAVIAR *fptr) : Constraint{fptr},
    domain{nullptr},
    dt{-1.0}, error_tolerance{1e-6},
    initialized{false}
{
  FC_OBJECT_INITIALIZE_INFO
}

M_shake::~M_shake () {}

bool M_shake::read (caviar::interpreter::Parser *parser) {
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
    } 
    else FC_ERR_UNDEFINED_VAR(t)

  }
  return in_file;
}

void M_shake::verify_settings () {
  FC_NULLPTR_CHECK(atom_data)
  atom_data->record_owned_position_old = true;

  FC_NULLPTR_CHECK(domain)
  if (dt <= 0.0) error->all (FC_FILE_LINE_FUNC, "dt have to be a positive number");      

}

void M_shake::step_part_I (int) {

}

void M_shake::step_part_II (int) {

}

void M_shake::step_part_III (int) {

  FC_OBJECT_VERIFY_SETTINGS

  bond_fix(); 


  // velocity_fix part
  // this fix has to be done only on the M-Shake molecules. If not, the normal
  // leap-frog step has to be enough.
  auto &vel = atom_data -> owned.velocity;
  auto &pos = atom_data -> owned.position;
  auto &pos_old = atom_data -> owned.position_old;
  auto &atomic_bond_index_vector = atom_data -> owned.atomic_bond_index_vector;
  for (unsigned int i=0; i<atomic_bond_index_vector.size(); i++) { 
    for (unsigned int j=0; j<atomic_bond_index_vector[i].size(); j++) { // XXX P.II
      const auto k = atomic_bond_index_vector[i][j];
      vel[k] = domain -> fix_distance(pos[k] - pos_old[k]) / dt ;			
    }
  }

}


void M_shake::bond_fix () {


  auto &pos = atom_data -> owned.position;
  auto &pos_old = atom_data -> owned.position_old;
  auto &atomic_bond_index_vector = atom_data -> owned.atomic_bond_index_vector;
  auto &atomic_bond_vector = atom_data -> owned.atomic_bond_vector;

  for (unsigned int i=0; i<atomic_bond_index_vector.size(); i++) { 

    auto Nc = atomic_bond_index_vector[i].size();
    if (Nc==0) continue;
	  std::vector<std::vector<double>> A (Nc, std::vector<double> (Nc,0));
    std::vector<std::vector<double>> M (Nc, std::vector<double> (Nc,0));  //inverse of matrix A


	  std::vector<double> C(Nc,0);

		double sum_err = 1.0;


    while(sum_err>error_tolerance) {

      for (unsigned int j=0; j<atomic_bond_vector[i].size(); j++) { 

        auto k1 = atomic_bond_vector[i][j].index_1, k2 = atomic_bond_vector[i][j].index_2;

        auto d = atomic_bond_vector[i][j].length;

        auto mass_inv_k1 = atom_data->owned.mass_inv[atom_data->owned.type[k1]];
        auto mass_inv_k2 = atom_data->owned.mass_inv[atom_data->owned.type[k2]];

        auto dr = domain -> fix_distance (pos[k1] - pos[k2]);

        auto r2 = dr * dr;

				C[j] = (r2 - d*d)/(4*dt*dt); 


        for (unsigned int jp=0; jp<atomic_bond_vector[i].size(); jp++) { 

          auto kp1 = atomic_bond_vector[i][jp].index_1, kp2 = atomic_bond_vector[i][jp].index_2;

          auto dr_old = domain -> fix_distance (pos[kp1] - pos[kp2]);

          auto r_r_old_dot = dr * dr_old;

				
          A[j][jp] = (mass_inv_k1*(delta(k1,kp1) - delta(k1,kp2)) 
                     +mass_inv_k2*(delta(k2,kp2) - delta(k2,kp1)) ) * r_r_old_dot;
					
				}
			}	



      int mi_err = caviar::matrix_inverse(A, M);

      if (mi_err==1)  error->all(FC_FILE_LINE_FUNC,"det==0");
      if (mi_err!=0)  error->all(FC_FILE_LINE_FUNC,"unknown error in matrix inverse");


      //lagrange multipliers
      std::vector<double> l(Nc, 0);		

			
      for(unsigned int t1=0;t1<Nc;t1++)
        for(unsigned int t2=0;t2<Nc;t2++)
          l[t1] += M[t1][t2]*C[t2];
      
	
	
      for (unsigned int j=0; j<atomic_bond_vector[i].size(); j++) { 
        int k1 = atomic_bond_vector[i][j].index_1, k2 = atomic_bond_vector[i][j].index_2;

        double mass_inv_k1 = atom_data->owned.mass_inv[atom_data->owned.type[k1]];
        double mass_inv_k2 = atom_data->owned.mass_inv[atom_data->owned.type[k2]];
	
				
        auto dr_old = domain -> fix_distance (pos_old[k1] - pos_old[k2]); 

        auto f_coef = -2.0 * dt * dt * l[j];

        auto fc = - f_coef * dr_old;

        pos[k1] -= fc * mass_inv_k1;
        pos[k2] += fc * mass_inv_k2;
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


