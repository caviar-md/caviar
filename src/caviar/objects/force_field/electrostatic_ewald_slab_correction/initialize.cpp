
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
#include "caviar/objects/domain.h"
#include "caviar/utility/macro_constants.h"

#include <cmath>
#include <iomanip>
#include <algorithm>

namespace caviar {
namespace objects {
namespace force_field {

void Electrostatic_ewald_slab_correction::initialize () {
  if (!calculated_once) {
    calculated_once = true;

    //const auto bc = domain->boundary_condition;
    //const auto bc_sum = bc.x + bc.y + bc.z;
    //if (bc_sum != 2)
    //  error->all(FC_FILE_LINE_FUNC, " only two "
    //            "dimensional periodic boundary condition is allowed.");

    const auto glx = domain->upper_global.x - domain->lower_global.x;
    const auto gly = domain->upper_global.y - domain->lower_global.y;
    const auto glz = domain->upper_global.z - domain->lower_global.z;

    if (slab_normal_axis == 0) {
      lx = gly; ly = glz; hz = glx;      
    } else if (slab_normal_axis == 1) {
      lx = glx; ly = glz; hz = gly;
    } else if (slab_normal_axis == 2) {
      lx = glx; ly = gly; hz = glz;
    }

    lx_ly_inv = 1.0/(lx*ly);

    if (hz > lz)
      error->all(FC_FILE_LINE_FUNC, "expected 'hz > lz' for slab geometry");
    slab_sum_e_coef = - FC_4PI * lx_ly_inv;
    make_slab_k_vectors();

    double l_xyz_inv = lx_ly_inv / lz;
    dipole_coef = k_electrostatic * FC_4PI * l_xyz_inv ;

  }

  calculate_dipole_sum();
  make_slab_chi_vectors();
}

void Electrostatic_ewald_slab_correction::calculate_dipole_sum() {

  const auto &pos = atom_data -> owned.position;
  auto ds = Vector<double> {0, 0, 0};
  for (unsigned int j=0;j<pos.size();++j) {
    const auto type_j = atom_data -> owned.type [j] ;
    const auto charge_j = atom_data -> owned.charge [ type_j ];   
    ds += charge_j* pos[j];
  }
  const auto dfv = - dipole_sum * dipole_coef ; // minus due to {-grad E}

  if (slab_normal_axis==0) {
    dipole_field_vector = Vector<double>{dfv.x, 0, 0};
    dipole_sum = Vector<double>{ds.x, 0, 0};
  }
  if (slab_normal_axis==1) {
    dipole_field_vector = Vector<double>{0, dfv.y, 0};
    dipole_sum = Vector<double>{0, ds.y, 0};
  }
  if (slab_normal_axis==2) {
    dipole_field_vector = Vector<double>{0, 0, dfv.z};
    dipole_sum = Vector<double>{0, 0, ds.z};
  }
}

void Electrostatic_ewald_slab_correction::make_slab_k_vectors () {

  kx.resize(kx_max);
  ky.resize(ky_max);
  kp.resize(kx_max*ky_max);

  kx_coef.resize(kx_max);
  ky_coef.resize(ky_max);
  kp_coef.resize(kx_max*ky_max);

  const auto lx_inv = 1.0/lx;
  const auto ly_inv = 1.0/ly;

#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for
#endif   
  for (auto ix = 0; ix < kx_max; ++ix) {
    kx[ix] = FC_2PI*(ix+1)*lx_inv;
    kx_coef[ix] = 1.0/(kx[ix]*(std::exp(kx[ix]*lz) - 1.0));
  }

#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for
#endif   
  for (auto iy = 0; iy < ky_max; ++iy) {
    ky[iy] = FC_2PI*(iy+1)*ly_inv;
    ky_coef[iy] = 1.0/(ky[iy]*(std::exp(ky[iy]*lz) - 1.0));
  }

  // XXX no parallel with openMP due  to external counter
  int ip = 0;
  for (auto ix = 0; ix <kx_max; ++ix) {
  for (auto iy = 0; iy <ky_max; ++iy) {
    kp[ip] = std::sqrt ((kx[ix]*kx[ix]) + (ky[iy]*ky[iy]));
    kp_coef[ip] = 1.0/(kp[ip]*(std::exp(kp[ip]*lz) - 1.0));
    ++ip;
  }
  }

}

void Electrostatic_ewald_slab_correction::make_slab_chi_vectors () {

///* XXX s1: To be checked which one is faster
  const std::vector<double> vec1d_0 {0.0, 0.0};
  const std::vector<std::vector<double>> vec2d_0 { vec1d_0, vec1d_0 };
  const std::vector<std::vector<std::vector<double>>> vec3d_0 { vec2d_0, vec2d_0 };
  chi_x.resize(kx_max, vec2d_0);
  chi_y.resize(ky_max, vec2d_0);
  chi_p.resize(ky_max*kx_max, vec3d_0);
  std::fill(chi_x.begin(), chi_x.end(), vec2d_0 );
  std::fill(chi_y.begin(), chi_y.end(), vec2d_0 );
  std::fill(chi_p.begin(), chi_p.end(), vec3d_0 );

//*/

/* // XXX s2: To be checked which one is faster
  chi_x.resize(kx_max);
  chi_y.resize(ky_max);
  chi_p.resize(ky_max*kx_max);

  for (int i = 0; i < kx_max; ++i) {
    std::cout << "i : " << i << std::endl;
    chi_x[i].resize(2);
  }
   std::cout << "here : "  << std::endl;
  for (int i = 0; i < kx_max; ++i)
  for (int j = 0; j < 2; ++j)
    chi_x[i][j].resize(2);

  for (auto && i :chi_y) i.resize(2);
  for (auto && i :chi_y) for (auto && j :i) j.resize(2);

  for (auto && i :chi_p) i.resize(2);
  for (auto && i :chi_p) for (auto && j :i) j.resize(2);
  for (auto && i :chi_p) for (auto && j :i)  for (auto && k :j) k.resize(2);

  for (auto && i :chi_x) for (auto && j :i) for (auto && k :j) k = 0.0;
  for (auto && i :chi_y) for (auto && j :i) for (auto && k :j) k = 0.0;
  for (auto && i :chi_p) for (auto && j :i) for (auto && k :j) for (auto && l :k) l = 0.0;
*/

  const auto &pos = atom_data->owned.position;
  const auto &type = atom_data->owned.type;
  const auto &charge = atom_data->owned.charge;
  const auto pos_size = pos.size();

  // XXX No openMP Parallel due to array reduction and boolean flag
  for (unsigned int j=0; j<pos_size; ++j) {
    const auto q = charge[ type[j] ];
    const auto p = give_slab_local_coordinates(pos[j]);


    bool do_iy_loop_once = true;
    int ip = 0;
    for (auto ix = 0; ix <kx_max; ++ix) {

      const double sinhkxz = std::sinh(kx[ix]*p.z);
      const double coshkxz = std::cosh(kx[ix]*p.z);

      const double sinkxx = std::sin(kx[ix]*p.x);
      const double coskxx = std::cos(kx[ix]*p.x);

      chi_x[ix][0][0] += q * sinhkxz * sinkxx;
      chi_x[ix][0][1] += q * sinhkxz * coskxx;
      chi_x[ix][1][0] += q * coshkxz * sinkxx;
      chi_x[ix][1][1] += q * coshkxz * coskxx;

      for (auto iy = 0; iy <ky_max; ++iy) {

        const double sinhkpz = std::sinh(kp[ip]*p.z);
        const double coshkpz = std::cosh(kp[ip]*p.z);


        const double sinkyy = std::sin(ky[iy]*p.y);
        const double coskyy = std::cos(ky[iy]*p.y);

        chi_p[ip][0][0][0] += q * sinhkpz * sinkxx * sinkyy;
        chi_p[ip][0][0][1] += q * sinhkpz * sinkxx * coskyy;
        chi_p[ip][0][1][0] += q * sinhkpz * coskxx * sinkyy;
        chi_p[ip][1][0][0] += q * coshkpz * sinkxx * sinkyy;
        chi_p[ip][0][1][1] += q * sinhkpz * coskxx * coskyy;
        chi_p[ip][1][0][1] += q * coshkpz * sinkxx * coskyy;
        chi_p[ip][1][1][0] += q * coshkpz * coskxx * sinkyy;
        chi_p[ip][1][1][1] += q * coshkpz * coskxx * coskyy;


        if (do_iy_loop_once) {
          const double sinhkyz = std::sinh(ky[iy]*p.z);
          const double coshkyz = std::cosh(ky[iy]*p.z);

          chi_y[iy][0][0] += q * sinhkyz * sinkyy;
          chi_y[iy][0][1] += q * sinhkyz * coskyy;
          chi_y[iy][1][0] += q * coshkyz * sinkyy;
          chi_y[iy][1][1] += q * coshkyz * coskyy;
        }

        ++ip;
      }

      do_iy_loop_once = false;
    }

  }


}

} //force_field
} //objects
} // namespace caviar

