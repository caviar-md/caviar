
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

namespace caviar {
namespace objects {
namespace force_field {

Vector<double> Electrostatic_ewald_slab_correction::field (const Vector<double> &r) {
  Vector<double> f_local{0,0,0};

// XXX Working Scheme of order N
// /*
// ======= kp sum
  const auto p = give_slab_local_coordinates(r); 
  int ip = 0;
  bool do_iy_loop_once = true;

// XXX no OpenMP parallel yet (due to boolean flag)
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:f_local)
#endif    
  for (auto ix = 0; ix <kx_max; ++ix) {


    const double sinhkxz = std::sinh(kx[ix]*p.z);
    const double coshkxz = std::cosh(kx[ix]*p.z);

    const double sinkxx = std::sin(kx[ix]*p.x);
    const double coskxx = std::cos(kx[ix]*p.x);

    const double sx_chi = - kx[ix]*( 
        - (+chi_x[ix][0][0] * coskxx * sinhkxz)
        - (-chi_x[ix][0][1] * sinkxx * sinhkxz)                            
        + (+chi_x[ix][1][0] * coskxx * coshkxz) 
        + (-chi_x[ix][1][1] * sinkxx * coshkxz) );

    const double sz_chi = - kx[ix]*( 
        - (+chi_x[ix][0][0] * coshkxz * sinkxx)
        - (+chi_x[ix][0][1] * coshkxz * coskxx)                            
        + (+chi_x[ix][1][0] * sinhkxz * sinkxx) 
        + (+chi_x[ix][1][1] * sinhkxz * coskxx) );

    f_local.x +=  kx_coef[ix] * sx_chi;
    f_local.z +=  kx_coef[ix] * sz_chi;


    for (auto iy = 0; iy <ky_max; ++iy) {

      const double sinhkpz = std::sinh(kp[ip]*p.z);
      const double coshkpz = std::cosh(kp[ip]*p.z);

      const double sinkyy = std::sin(ky[iy]*p.y);
      const double coskyy = std::cos(ky[iy]*p.y);


      const double sx_chi = - kx[ix]*( 
          - (+chi_p[ip][0][0][0] * coskxx * sinhkpz * sinkyy) 
          - (+chi_p[ip][0][0][1] * coskxx * sinhkpz * coskyy) 
          - (-chi_p[ip][0][1][0] * sinkxx * sinhkpz * sinkyy) 
          + (+chi_p[ip][1][0][0] * coskxx * coshkpz * sinkyy) 
          - (-chi_p[ip][0][1][1] * sinkxx * sinhkpz * coskyy) 
          + (+chi_p[ip][1][0][1] * coskxx * coshkpz * coskyy) 
          + (-chi_p[ip][1][1][0] * sinkxx * coshkpz * sinkyy) 
          + (-chi_p[ip][1][1][1] * sinkxx * coshkpz * coskyy) );

      const double sy_chi = - ky[iy]*( 
          - (+chi_p[ip][0][0][0] * coskyy * sinhkpz * sinkxx) 
          - (-chi_p[ip][0][0][1] * sinkyy * sinhkpz * sinkxx) 
          - (+chi_p[ip][0][1][0] * coskyy * sinhkpz * coskxx) 
          + (+chi_p[ip][1][0][0] * coskyy * coshkpz * sinkxx) 
          - (-chi_p[ip][0][1][1] * sinkyy * sinhkpz * coskxx) 
          + (-chi_p[ip][1][0][1] * sinkyy * coshkpz * sinkxx) 
          + (+chi_p[ip][1][1][0] * coskyy * coshkpz * coskxx) 
          + (-chi_p[ip][1][1][1] * sinkyy * coshkpz * coskxx) );

      const double sz_chi = - kp[ip]*(
          - (+chi_p[ip][0][0][0] * coshkpz * sinkxx* sinkyy) 
          - (+chi_p[ip][0][0][1] * coshkpz * sinkxx* coskyy) 
          - (+chi_p[ip][0][1][0] * coshkpz * coskxx* sinkyy) 
          + (+chi_p[ip][1][0][0] * sinhkpz * sinkxx* sinkyy) 
          - (+chi_p[ip][0][1][1] * coshkpz * coskxx* coskyy) 
          + (+chi_p[ip][1][0][1] * sinhkpz * sinkxx* coskyy) 
          + (+chi_p[ip][1][1][0] * sinhkpz * coskxx* sinkyy) 
          + (+chi_p[ip][1][1][1] * sinhkpz * coskxx* coskyy) );

      f_local.x +=  2.0 * kp_coef[ip] * sx_chi;
      f_local.y +=  2.0 * kp_coef[ip] * sy_chi;
      f_local.z +=  2.0 * kp_coef[ip] * sz_chi;

      if (do_iy_loop_once) {
        const double sinhkyz = std::sinh(ky[iy]*p.z);
        const double coshkyz = std::cosh(ky[iy]*p.z);
    
        const double syy_chi = - ky[iy]*( 
            - (+chi_y[iy][0][0] * coskyy * sinhkyz) 
            - (-chi_y[iy][0][1] * sinkyy * sinhkyz) 
            + (+chi_y[iy][1][0] * coskyy * coshkyz) 
            + (-chi_y[iy][1][1] * sinkyy * coshkyz) );

        const double syz_chi = - ky[iy]*( 
            - (+chi_y[iy][0][0] * sinkyy * coshkyz) 
            - (+chi_y[iy][0][1] * coskyy * coshkyz) 
            + (+chi_y[iy][1][0] * sinkyy * sinhkyz) 
            + (+chi_y[iy][1][1] * coskyy * sinhkyz) );

        f_local.y +=  ky_coef[iy] * syy_chi;
        f_local.z +=  ky_coef[iy] * syz_chi;
      }
      ++ip;
    }      
    do_iy_loop_once = false;
  }

  f_local *= 2.0 * slab_sum_e_coef * k_electrostatic;

// */
//


// XXX Working Scheme of order N^2
 /*
  const auto &type = atom_data->owned.type;
  const auto &charge = atom_data->owned.charge;
  const auto &pos = atom_data->owned.position;
  const auto pos_size = pos.size();  

  for (auto j = 0; j < pos_size; ++j) {
    const auto pos_j = pos[j];
    const auto pos_ij = r - pos_j; // XXX ? or its negative ?
    const auto plij = give_slab_local_coordinates (pos_ij);
    const auto charge_j = charge[type[j]];
    


// ======= kp sum

    int ip = 0;
    for (auto ix = 0; ix <kx_max; ++ix) {
    for (auto iy = 0; iy <ky_max; ++iy) {

    const double sinhkpz = std::sinh(kp[ip]*plij.z);
    const double coshkpz = std::cosh(kp[ip]*plij.z);

    const double sinkxx = std::sin(kx[ix]*plij.x);
    const double coskxx = std::cos(kx[ix]*plij.x);

    const double sinkyy = std::sin(ky[iy]*plij.y);
    const double coskyy = std::cos(ky[iy]*plij.y);


    f_local.x +=  2.0 * charge_j * kp_coef[ip] * kx[ix] * coshkpz * sinkxx * coskyy;
    f_local.y +=  2.0 * charge_j * kp_coef[ip] * ky[iy] * coshkpz * coskxx * sinkyy;
    f_local.z += -2.0 * charge_j * kp_coef[ip] * kp[ip] * sinhkpz * coskxx * coskyy;
    ++ip;
  }
  }

//========== kx sum

  for (auto ix=0; ix < kx_max; ++ix) {

    const double sinhkxz = std::sinh(kx[ix]*plij.z);
    const double coshkxz = std::cosh(kx[ix]*plij.z);

    const double sinkxx = std::sin(kx[ix]*plij.x);
    const double coskxx = std::cos(kx[ix]*plij.x);


    f_local.x +=  charge_j * kx_coef[ix] * kx[ix] * coshkxz * sinkxx;
    f_local.z += -charge_j * kx_coef[ix] * kx[ix] * sinhkxz * coskxx;
  }

//===========ky sum

  for (auto iy=0; iy < ky_max; ++iy) {

    const double sinhkyz = std::sinh(ky[iy]*plij.z);
    const double coshkyz = std::cosh(ky[iy]*plij.z);

    const double sinkyy = std::sin(ky[iy]*plij.y);
    const double coskyy = std::cos(ky[iy]*plij.y);


    f_local.y +=  charge_j * ky_coef[iy] *  ky[iy] * coshkyz * sinkyy;
    f_local.z += -charge_j * ky_coef[iy] *  ky[iy] * sinhkyz * coskyy;
  }

  }

  f_local *= 2.0 * slab_sum_e_coef * k_electrostatic;
// */
//

  return give_slab_global_coordinates(f_local) + dipole_field_vector;
}

Vector<double> Electrostatic_ewald_slab_correction::field (int i) {
  return field (atom_data->owned.position[i]);
}

Vector<double> Electrostatic_ewald_slab_correction::dipole_field () {
  return dipole_field_vector;
}

} //force_field
} //objects
} // namespace caviar

