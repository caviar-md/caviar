
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

#include "caviar/objects/force_field/electrostatic_ewald_k.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/macro_constants.h"

#include <complex>
#include <cmath>
#include <iomanip>

namespace caviar {

namespace force_field {


double Electrostatic_ewald_k::potential (const Vector<double> &r) {
  double p = k_space_potential(r);
  if (dipole)
    p += dipole_potential(r);
  return p;
}

double Electrostatic_ewald_k::potential (const int i) {
  double p = k_space_potential(i);
  if (dipole)
    p += dipole_potential(i);
  return p;
}

double Electrostatic_ewald_k::k_space_potential (const Vector<double> &r) {
// XXX Working scheme
/*
  double sum_k = 0;

  static std::complex<double> ii(0.0, 1.0);    
  const auto &pos = atom_data -> owned.position;

  for (unsigned int j=0;j<pos.size();++j) {
    const auto type_j = atom_data -> owned.type [j] ;
    const auto charge_j = atom_data -> owned.charge [ type_j ];
    const auto r_ij = r - pos[j];
    std::complex<double> rho (0,0);
    for (int k = 0; k < n_k_vectors; ++k) {
      rho +=  field_k_coef[k] * std::exp(-ii*(k_vector[k]*r_ij));
    }

    const double rho_norm = std::real(rho);
    sum_k += charge_j * rho_norm;

  }

  return FC_4PI * l_xyz_inv * sum_k * k_electrostatic;
*/

//Re(S_j S_k e(-ik(ri-rj))) = Re(S_k S_j e(-ik(ri-rj))) = Re(S_k e(-ik ri) S_j e(+ik(rj)))

// XXX Working scheme with seperation of 'r' and 'pos_j' vectors
/*
  double sum_k = 0;

  static std::complex<double> ii(0.0, 1.0);    
  const auto &pos = atom_data -> owned.position;

  for (int k = 0; k < n_k_vectors; ++k) {
    const auto k_vector_k = k_vector[k];
    const auto field_k_coef_k = field_k_coef[k];

    std::complex<double> rho (0,0);

    for (unsigned int j=0;j<pos.size();++j) { 
      const auto type_j = atom_data -> owned.type [j] ;
      const auto charge_j = atom_data -> owned.charge [ type_j ];

      rho +=  charge_j * std::exp(ii*(k_vector_k*pos[j]));

    }

    rho *= std::exp(-ii*(k_vector_k*r));

    const double rho_norm = field_k_coef_k*std::real(rho);

    sum_k += rho_norm;
  }
  return FC_4PI * l_xyz_inv * sum_k * k_electrostatic;
*/

// XXX Fastest Working scheme
// /*
  double sum_k = 0;

  static std::complex<double> ii(0.0, 1.0);    

#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:sum_k)
#endif        
  for (int k = 0; k < n_k_vectors; ++k) {
    const auto k_vector_k = k_vector[k];
    const auto field_k_coef_k = field_k_coef[k];

    auto rho = potential_k_coef_cmplx[k];

    rho *= std::exp(-ii*(k_vector_k*r));

    const double rho_norm = field_k_coef_k*std::real(rho);

    sum_k += rho_norm;
  }
  return FC_4PI * l_xyz_inv * sum_k * k_electrostatic; // XXX add dipole potential here
// */

}

double Electrostatic_ewald_k::k_space_potential (const int i) {
  return k_space_potential(atom_data-> owned.position[i]);
}

double Electrostatic_ewald_k::dipole_potential (const Vector<double> &r) {
  return 0*r.x; // XXX
}

double Electrostatic_ewald_k::dipole_potential (const int i) {
  return 0*i; // XXX is it correct? Since the dipole field is position independent
            // We also expect that the potential be periodic. Do we?
}


} //force_field

} // namespace caviar

