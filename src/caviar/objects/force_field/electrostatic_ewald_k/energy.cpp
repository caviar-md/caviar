
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

#include <cmath>
#include <iomanip>
#include <complex>

namespace caviar {

namespace force_field {


double Electrostatic_ewald_k::self_energy () {
  const auto &pos = atom_data -> owned.position;    
  double sum_j {0};
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:sum_j)
#endif  
  for (unsigned int j=0;j<pos.size();++j) {
    const auto type_j = atom_data -> owned.type [j] ;  
    const auto charge_j = atom_data -> owned.charge [ type_j ];
    sum_j += charge_j* charge_j;  //
  }
  return - FC_PIS_INV * alpha*sum_j* k_electrostatic;
}

double Electrostatic_ewald_k::dipole_energy () {
//XXX paper's Energy folmula
  return 0.5*dipole_coef*(dipole_sum*dipole_sum);
// XXX another idea : Dipole_Energy = -Dipole_field*Dipole_value
// There's a 0.5 factor needed to get paper's formula
//  return -0.5*(dipole_field_vector*dipole_sum);
}

double Electrostatic_ewald_k::k_space_energy () {
// Working Scheme using energy formula. Of order 'N'
/*
  const auto &pos = atom_data -> owned.position;    
  const std::complex<double> ii(0.0, 1.0);    
  double sum_k = 0.0;

  for (int i = 0; i<n_k_vectors; ++i) {

    std::complex<double> rho (0,0);
    for (unsigned int j=0;j<pos.size();++j) {
      const auto type_j = atom_data -> owned.type [j] ;
      const auto charge_j = atom_data -> owned.charge [ type_j ];

      rho += charge_j * std::exp(-ii*(k_vector[i]*pos[j]));
    }
    const double rho_norm = std::abs(rho * std::conj(rho));
    sum_k += field_k_coef[i] * rho_norm;
  }    

  return 0.5 * l_xyz_inv * FC_4PI * sum_k * k_electrostatic;
*/


// XXX Working scheme using potential function. Order depends on the order
// of the potential. It can be changing from N^2 to N..
///*
  const auto &pos = atom_data -> owned.position;    
  double energy_k = 0;
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:energy_k)
#endif  
  for (unsigned int j=0;j<pos.size();++j) {
    const auto type_j = atom_data -> owned.type [j] ;  
    const auto charge_j = atom_data -> owned.charge [ type_j ];
    energy_k += charge_j * k_space_potential(pos[j]); // working
//    energy_k += charge_j * k_space_potential(j); // working
  }
  return 0.5 * energy_k;   
//*/

}

double Electrostatic_ewald_k::energy () {

  double e = k_space_energy() + self_energy();
  if (dipole) 
    e += dipole_energy();
  return e;
}

} //force_field

} // namespace caviar

