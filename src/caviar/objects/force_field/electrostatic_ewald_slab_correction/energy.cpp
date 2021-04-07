
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

#include "caviar/objects/force_field/electrostatic_ewald_slab_correction.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/macro_constants.h"

#include <cmath>
#include <iomanip>
#include <complex>

#define FC_POW2(X) (X*X)

namespace caviar {
namespace objects {
namespace force_field {


double Electrostatic_ewald_slab_correction::energy () {

  double sum_e = 0.0;

// XXX Optimized working Scheme
// /* 

// ======= kp sum

#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:sum_e)
#endif     
  for (auto ip=0; ip < kx_max*ky_max; ++ip) {

    const double sum_chi = 
        - FC_POW2(chi_p[ip][0][0][0]) - FC_POW2(chi_p[ip][0][0][1]) 
        - FC_POW2(chi_p[ip][0][1][0]) + FC_POW2(chi_p[ip][1][0][0]) 
        - FC_POW2(chi_p[ip][0][1][1]) + FC_POW2(chi_p[ip][1][0][1]) 
        + FC_POW2(chi_p[ip][1][1][0]) + FC_POW2(chi_p[ip][1][1][1]);


   sum_e +=  2.0 * kp_coef[ip] * sum_chi;
  }

//========== kx sum
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:sum_e)
#endif     
  for (auto ix=0; ix < kx_max; ++ix) {

    const double sum_chi = 
        - FC_POW2(chi_x[ix][0][0]) 
        - FC_POW2(chi_x[ix][0][1]) 
        + FC_POW2(chi_x[ix][1][0]) 
        + FC_POW2(chi_x[ix][1][1]);


    sum_e +=  kx_coef[ix] * sum_chi;
  }

//===========ky sum
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:sum_e)
#endif     
  for (auto iy=0; iy < ky_max; ++iy) {

    const double sum_chi =
        - FC_POW2(chi_y[iy][0][0])
        - FC_POW2(chi_y[iy][0][1]) 
        + FC_POW2(chi_y[iy][1][0])
        + FC_POW2(chi_y[iy][1][1]);

    sum_e += ky_coef[iy] * sum_chi;
  }

  sum_e *= k_electrostatic * slab_sum_e_coef;
// */
//


// XXX working Scheme of Order N^2
 /*

  const auto &pos = atom_data->owned.position;
  const auto &type = atom_data->owned.type;
  const auto &charge = atom_data->owned.charge;
  const auto pos_size = pos.size();
  for (unsigned int i = 0; i < pos_size; ++i) {
  for (unsigned int j = 0; j < pos_size; ++j) {

    const auto p = give_slab_local_coordinates(pos[i]-pos[j]); 
    const double q_ij = charge[ type[i] ] * charge[ type[j] ];

    double sum_kp = 0;
    int ip = 0;
    for (int ix=0; ix < kx_max; ++ix) {
    for (int iy=0; iy < ky_max; ++iy) {
      sum_kp += kp_coef[ip] * std::cosh(kp[ip]*p.z)*std::cos(kx[ix]*p.x)*std::cos(ky[iy]*p.y);
      ++ip;
    }
    }
    sum_e += q_ij *  2.0 * sum_kp;


    double sum_kx = 0;
    for (int ix=0; ix < kx_max; ++ix) {
      sum_kx += kx_coef[ix] *std::cosh(kx[ix]*p.z)*std::cos(kx[ix]*p.x);
    }
    sum_e += q_ij * sum_kx;


    double sum_ky = 0;
    for (int iy=0; iy < ky_max; ++iy) {
      sum_ky += ky_coef[iy] *std::cosh(ky[iy]*p.z)*std::cos(ky[iy]*p.y);
    }
    sum_e += q_ij * sum_ky;

  }
  }
  sum_e *= k_electrostatic * slab_sum_e_coef;
 */

// XXX working Scheme using potential function
 /* 

  const auto &pos = atom_data->owned.position;
  const auto &type = atom_data->owned.type;
  const auto &charge = atom_data->owned.charge;
  const auto pos_size = pos.size();
  for (int i = 0; i < pos_size; ++i) {
    const auto charge_i = charge[ type [i] ];

    // XXX does it need a '0.5' coefficient ?
    sum_e +=  charge_i * potential(pos[i]);
  } 

 */

  return sum_e + dipole_energy();
}

double Electrostatic_ewald_slab_correction::dipole_energy () {
 //XXX paper's Energy folmula
  return 0.5 * dipole_coef * (dipole_sum*dipole_sum);
// XXX another idea : Dipole_Energy = -Dipole_field*Dipole_value
// There's a 0.5 factor needed to get paper's formula
//  return -0.5*(dipole_field_vector*dipole_sum);
}

} //force_field
} //objects
} // namespace caviar


