
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATIC_EWALD_K_H
#define CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATIC_EWALD_K_H

#include "caviar/objects/force_field.h"

#include <complex>

namespace caviar {
namespace objects {
 
namespace force_field {

/**
 * This class is electrostatic ewald in k-space force-field.
 * 
 * 
 */
class Electrostatic_ewald_k : public Force_field {
public:
  Electrostatic_ewald_k (class CAVIAR *);
  ~Electrostatic_ewald_k () {};
  
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
  double energy ();

public:

  double self_energy ();
  double dipole_energy ();
  double k_space_energy ();

  double potential (const Vector<double> &);
  double potential (const int);
  double k_space_potential (const Vector<double> &);
  double k_space_potential (const int);
  double dipole_potential (const Vector<double> &);
  double dipole_potential (const int);

  Vector<double> field (const Vector<double> &);
  Vector<double> field (const int);
  Vector<double> k_space_field (const Vector<double> &);
  Vector<double> k_space_field (const int);
  Vector<double> dipole_field ();


  int n_k_vectors;
  bool initialized, calculated_once;
  void initialize();
  void make_k_vectors();
  void calculate_alpha();
  void calculate_dipole_sum();

  double k_electrostatic, alpha;
  Vector<double> external_field;

// simulation box lengths and its product
  double lx,lx_inv, ly,ly_inv, lz,lz_inv, l_xyz_inv;

  int kx_max, ky_max, kz_max;
//  std::vector<std::vector<std::vector<double>>> k_coef;

  std::vector<Vector<double>> k_vector;
  std::vector<double> k_vector_sq;
  std::vector<double> field_k_coef;//, potential_k_coef;

  std::vector<std::complex<double>> potential_k_coef_cmplx; // depends on the positions...
  void calculate_potential_k_coef_cmplx(); // ...of the particles in each time step

  bool dipole;
  double epsilon_dipole;
  double dipole_coef; // defined as k_electrostatic * 4PI/(1+2e')L^3
  Vector<double> dipole_field_vector, dipole_sum;  // dipole_sum: Sum_j(q_j vec(r_j))
  
  bool slab_geometry;
  int slab_normal_axis;
  double slab_lz;


};

} //force_field
} //objects
} // namespace caviar

#endif
