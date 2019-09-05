
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
#include "caviar/objects/domain.h"
#include "caviar/utility/macro_constants.h"

#include <cmath>
#include <iomanip>

namespace caviar {
namespace objects {
namespace force_field {

void Electrostatic_ewald_k::initialize () {
  if (!calculated_once) {
    calculated_once = true;

    const auto bc = domain->boundary_condition;
    const auto bc_sum = bc.x + bc.y + bc.z;
    if (bc_sum < 2) output->warning("Electrostatic_ewald_k:: domain boundary "
        "condition has to be periodic in at least two dimensions. It doesn't "
        "affect 'Electrostatic_ewald_k' but it will affect 'Electrostatic_ewald_r'.");
    if (bc_sum == 2 && slab_geometry == false) {
      output->warning(" Note that an ewald slab correction is needed to gain " 
                      "correct results when there's two dim. periodic boundary condition.");
    }

    lx = domain->upper_global.x - domain->lower_global.x;
    ly = domain->upper_global.y - domain->lower_global.y;
    lz = domain->upper_global.z - domain->lower_global.z;

    if (slab_geometry) {
      if (slab_normal_axis == 0) lx = slab_lz;
      if (slab_normal_axis == 1) ly = slab_lz;
      if (slab_normal_axis == 2) lz = slab_lz;

      bool inconsistency = true;
      if ((bc.x==0) && (slab_normal_axis==0)) inconsistency = false;
      if ((bc.y==0) && (slab_normal_axis==1)) inconsistency = false;
      if ((bc.z==0) && (slab_normal_axis==2)) inconsistency = false;

      if (inconsistency)
        output->warning(" There's an inconsistency between 'slab_normal_axis' " 
                        "and boundary condition. 'slab_normal_axis' should be non-periodic");
    }

    lx_inv = 1.0/lx;  ly_inv = 1.0/ly;  lz_inv = 1.0/lz;
    l_xyz_inv = 1.0 / (lx*ly*lz);

    make_k_vectors();

    if (dipole) {
      dipole_coef = k_electrostatic * FC_4PI * l_xyz_inv / (1+2.0*epsilon_dipole);
    }
  }

  if (dipole) 
    calculate_dipole_sum();

  calculate_potential_k_coef_cmplx();
}

void Electrostatic_ewald_k::calculate_dipole_sum() {
    const auto &pos = atom_data -> owned.position;
    dipole_sum = Vector<double> {0, 0, 0};
    for (unsigned int j=0;j<pos.size();++j) {
      const auto type_j = atom_data -> owned.type [j] ;
      const auto charge_j = atom_data -> owned.charge [ type_j ];   
      dipole_sum += charge_j* pos[j];
    }
    dipole_field_vector = - dipole_sum*dipole_coef ;
}

void Electrostatic_ewald_k::calculate_potential_k_coef_cmplx() {
  static std::complex<double> ii(0.0, 1.0);    
  const auto &pos = atom_data -> owned.position;
  const auto pos_size = pos.size();
  potential_k_coef_cmplx.resize(n_k_vectors);

  for (int k = 0; k < n_k_vectors; ++k) {
    const auto k_vector_k = k_vector[k];

    std::complex<double> rho (0,0);

    for (unsigned int j=0;j<pos_size;++j) { 
      const auto type_j = atom_data -> owned.type [j] ;
      const auto charge_j = atom_data -> owned.charge [ type_j ];

      rho +=  charge_j * std::exp(ii*(k_vector_k*pos[j]));

    }
    potential_k_coef_cmplx[k] = rho;
  }
}


void Electrostatic_ewald_k::calculate_alpha () {
/*
  output->info("Electrostatic_ewald_k: calculate_alpha");
  int sum_of_setup = static_cast<int>(by_alpha) + static_cast<int>(by_accuracy)
                    + static_cast<int>(by_cutoff_r) ;
  if (sum_of_setup > 1) error->all("Electrostatic_ewald_k::calculate: accuracy, cutoff and alpha cannot be set together.");
  if (sum_of_setup == 0) error->all("Electrostatic_ewald_k::calculate: accuracy or cutoff or alpha has to be set.");
  if (by_alpha) {
    error->all("Electrostatic_ewald_k::calculate_alpha: not developed yet. a  cutoff formula");
  }
  if (by_cutoff_r) {
    // W. Smith, Information Quarterly for Computer Simulation of Condensed Phases, 21 (1986) 37
    alpha = 3.5 / cutoff_r; 
    error->all("Electrostatic_ewald_k::calculate_alpha: not developed yet");
  }
  if (by_accuracy) {
    if (domain==nullptr) error->all("Electrostatic_ewald_k::calculate: domain_set = nullptr");
    if (neighborlist==nullptr) error->all("Electrostatic_ewald_k::calculate: neighborlist = nullptr");
    cutoff_r = neighborlist -> cutoff;
    const auto no_bins = neighborlist -> no_bins;

    double ln_accuracy = std::log(accuracy);
    double no_bins_mean = (1.0/3.0)*(no_bins.x + no_bins.y + no_bins.z);
    // J. Perram, H. Petersen and S. De Leeuw, Mol. Phys. 65 (1988) 875.
    // there's a modification by no basis, which is the mean of no_bins that used.
    // in the paper there's : alpha = no_bins * std::sqrt(...)
    alpha =  no_bins_mean * std::sqrt (-ln_accuracy); 
    double P_inv = 1.0 / M_PI;
    kx_max = - no_bins.x * P_inv * ln_accuracy;
    ky_max = - no_bins.y * P_inv * ln_accuracy;
    kz_max = - no_bins.z * P_inv * ln_accuracy;
    // 
  }*/
}

void Electrostatic_ewald_k::make_k_vectors () {

  k_vector.clear();
  k_vector_sq.clear();
  field_k_coef.clear();

  n_k_vectors = (2*kx_max)*(2*ky_max)*(2*kz_max);

  //k_vector.reserve(n_k_vectors);
  //k_vector_sq.reserve(n_k_vectors);
  //field_k_coef.reserve(n_k_vectors);

  const auto alpha_sq = alpha*alpha;
  const auto four_alpha2_inv = 1.0/(4.0*alpha_sq);

  for (auto kx = -kx_max; kx <=kx_max; ++kx) {
  for (auto ky = -ky_max; ky <=ky_max; ++ky) {
  for (auto kz = -kz_max; kz <=kz_max; ++kz) {

    Vector<double> k_vec {FC_2PI*kx*lx_inv, FC_2PI*ky*ly_inv, FC_2PI*kz*lz_inv};

    const auto k_vec_sq = k_vec * k_vec;
    if (k_vec_sq == 0) continue;

    auto k_vec_sq_inv = 1.0 / k_vec_sq;

    k_vector.emplace_back(k_vec);
    k_vector_sq.emplace_back(k_vec_sq);
    field_k_coef.emplace_back(k_vec_sq_inv * std::exp(-k_vec_sq*four_alpha2_inv));

  }}}

}

} //force_field
} //objects
} // namespace caviar

