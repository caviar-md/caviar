
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

#include "caviar/objects/force_field/electrostatic_p3m.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/macro_constants.h"

#include <cmath>
#include <iomanip>
#include <complex>

namespace caviar {
namespace objects {
namespace force_field {

Electrostatic_p3m::Electrostatic_p3m (CAVIAR *fptr) : Force_field{fptr}, 
k_electrostatic{1.0}
 {
  FC_OBJECT_INITIALIZE_INFO
  error->all (FC_FILE_LINE_FUNC, "P3M is NOT IMPLEMENTED YET.");      
  initialized = false;
  calculated_once = false;
  dipole = false;
  epsilon_dipole = 1.0;
  kx_max = 1;  ky_max = 1;  kz_max = 1;
  slab_geometry = false;
  slab_normal_axis = 2;
  slab_lz = 1.0;
}

bool Electrostatic_p3m::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"k_max")) {
      int k_max = 1;
      GET_OR_CHOOSE_A_INT(k_max,"","")
      if (k_max < 0) error->all (FC_FILE_LINE_FUNC_PARSE, "'k_max' have to non-negative.");      
      kx_max = ky_max = kz_max = k_max;
    } else if (string_cmp(t,"kx_max")) {
      GET_OR_CHOOSE_A_INT(kx_max,"","")
      if (kx_max < 0) error->all (FC_FILE_LINE_FUNC_PARSE, "'kx_max' have to non-negative.");      
    } else if (string_cmp(t,"ky_max")) {
      GET_OR_CHOOSE_A_INT(ky_max,"","")
      if (ky_max < 0) error->all (FC_FILE_LINE_FUNC_PARSE, "'ky_max' have to non-negative.");      
    } else if (string_cmp(t,"kz_max")) {
      GET_OR_CHOOSE_A_INT(kz_max,"","")
      if (kz_max < 0) error->all (FC_FILE_LINE_FUNC_PARSE, "'kz_max' have to non-negative.");      
    } else if (string_cmp(t,"slab_normal_axis")) {
      slab_geometry = true;
      const auto t = parser->get_val_token();
      if (t.kind == caviar::interpreter::Kind::identifier) {
        const auto ts = t.string_value;
        slab_normal_axis = -1;
        if (string_cmp(ts, "x")) slab_normal_axis = 0;
        if (string_cmp(ts, "y")) slab_normal_axis = 1;
        if (string_cmp(ts, "z")) slab_normal_axis = 2;
        if (slab_normal_axis == -1)
          error->all (FC_FILE_LINE_FUNC_PARSE, "expected 'x' or 'y' or 'z' for axis."); 
      }
    } else if (string_cmp(t,"k_electrostatic")) {
      GET_OR_CHOOSE_A_REAL(k_electrostatic,"","")    
      if (k_electrostatic < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "'k_electrostatic' has to be non-negative.");            
    } else if (string_cmp(t,"epsilon_dipole")) {
      GET_OR_CHOOSE_A_REAL(epsilon_dipole,"","")    
      if (epsilon_dipole < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "'epsilon_dipole' has to be non-negative.");            
      dipole = true;
    } else if (string_cmp(t,"slab_lz")||string_cmp(t,"slab_thickness")) {
      slab_geometry = true;
      GET_OR_CHOOSE_A_REAL(slab_lz,"","")
      if (slab_lz < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "slab_lz have to non-negative.");  
    } else if (string_cmp(t,"alpha")) {
      GET_OR_CHOOSE_A_REAL(alpha,"","")
      if (alpha < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "alpha have to non-negative.");  
    } else if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"set_domain") || string_cmp(t,"domain")) {
      FIND_OBJECT_BY_NAME(domain,it)
      domain = object_container->domain[it->second.index];
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  
  return in_file;
}

void Electrostatic_p3m::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)
  FC_NULLPTR_CHECK(domain)
}


void Electrostatic_p3m::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS

  initialize();

/* XXX working scheme but of order N^2
  const auto &pos = atom_data -> owned.position;    
  for (unsigned int i=0;i<pos.size();++i) {
  const auto type_i = atom_data -> owned.type [i] ;
  const auto charge_i = atom_data -> owned.charge [ type_i ];    
  const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];    
  Vector<double> sum_j {0,0,0};
  for (unsigned int j=0;j<pos.size();++j) {

    const auto type_j = atom_data -> owned.type [j] ;
    const auto charge_j = atom_data -> owned.charge [ type_j ];    
    const auto r_ij = pos[i] - pos[j];
    Vector<double> sum_k {0,0,0};
    for (int k = 0; k<n_k_vectors; ++k) {
      sum_k += FC_4PI * field_k_coef[k] * k_vector[k] * std::sin(k_vector[k]*r_ij);
    }
    sum_j += charge_j * sum_k;
  }
  const auto force =  k_electrostatic * charge_i * l_xyz_inv * sum_j;
  atom_data -> owned.acceleration[i] += force * mass_inv_i;
  }
*/

// XXX Working scheme of the order N. Best one.
///*
  const auto &pos = atom_data -> owned.position;    
  static std::complex<double> ii(0.0, 1.0);    


  for (int k = 0; k<n_k_vectors; ++k) {

    const auto sum_j_c = std::conj(potential_k_coef_cmplx[k]);
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
    for (unsigned int i=0;i<pos.size();++i) {
      const auto type_i = atom_data -> owned.type [i] ;
      const auto charge_i = atom_data -> owned.charge [ type_i ];    
      const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];    
  
      const double sum_j = std::imag(sum_j_c * std::exp(ii*(k_vector[k]*pos[i])));
      const auto sum_k = FC_4PI * field_k_coef[k] * k_vector[k] * sum_j ;

      auto force =  k_electrostatic * charge_i * l_xyz_inv * sum_k;

      if (dipole) {
        force += charge_i * dipole_field_vector;
      }

      atom_data -> owned.acceleration[i] += force * mass_inv_i;     

    }
  }
//*/

// XXX Working scheme of the order N.
/*
  const auto &pos = atom_data -> owned.position;    
  static std::complex<double> ii(0.0, 1.0);    

  for (unsigned int i=0;i<pos.size();++i) {
    const auto type_i = atom_data -> owned.type [i] ;
    const auto charge_i = atom_data -> owned.charge [ type_i ];    
    const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];    

    std::complex<double> sum_kx (0.0, 0.0);
    std::complex<double> sum_ky (0.0, 0.0);
    std::complex<double> sum_kz (0.0, 0.0);

    for (int k = 0; k<n_k_vectors; ++k) {

      const auto sum_j_c = std::conj(potential_k_coef_cmplx[k]);
      const auto c = field_k_coef[k] * sum_j_c * std::exp(ii*(k_vector[k]*pos[i]));   
      sum_kx += k_vector[k].x * c;
      sum_ky += k_vector[k].y * c;
      sum_kz += k_vector[k].z * c;

    }  

    const double sum_jx = std::imag(sum_kx);
    const double sum_jy = std::imag(sum_ky);
    const double sum_jz = std::imag(sum_kz);

    Vector<double> sv {sum_jx, sum_jy, sum_jz};

    const auto sum = FC_4PI *  sv ;

    const auto force =  k_electrostatic * charge_i * l_xyz_inv * sum;
    atom_data -> owned.acceleration[i] += force * mass_inv_i;

  }
 
*/

// XXX Scheme using field functions.
 /*
  const auto &pos = atom_data -> owned.position;
  const unsigned pos_size = pos.size();

  for (unsigned i = 0; i < pos_size; ++i) {
    const auto pos_i = atom_data->owned.position [i];
    const auto type_i = atom_data -> owned.type [ i ];      
    const auto charge_i = atom_data -> owned.charge [ type_i ];      
    const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];      

    const auto force = charge_i * field (pos_i);  // working
//    const auto force = charge_i * field (i);  // working
    atom_data -> owned.acceleration[i] += force * mass_inv_i;

  }

 */

}

} //force_field
} //objects
} // namespace caviar

