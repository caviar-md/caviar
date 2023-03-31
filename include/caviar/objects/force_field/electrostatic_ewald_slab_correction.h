
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATICEWALDSLABCORRECTION_H
#define CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATICEWALDSLABCORRECTION_H

#include "caviar/objects/force_field.h"

#include <complex>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  /**
   * This class has electric double layer correction ELC for slab geometry.
   * It should be used alongside ewald_k and ewald_r methods
   *
   */
  class Electrostatic_ewald_slab_correction : public Force_field
  {
  public:
    Electrostatic_ewald_slab_correction(class CAVIAR *);
    ~Electrostatic_ewald_slab_correction(){};

    bool read(class caviar::interpreter::Parser *);
    void verify_settings();
    void calculate_acceleration();

    double energy();

  public:
    int slab_normal_axis;
    double potential(const Vector<double> &);
    double potential(const int);

    Vector<double> field(const Vector<double> &);
    Vector<double> field(const int);

    inline Vector<double> give_slab_local_coordinates(const Vector<double> &vg)
    {
      Vector<double> vl;
      if (slab_normal_axis == 0)
      {
        vl.x = vg.y;
        vl.y = vg.z;
        vl.z = vg.x;
      }
      else if (slab_normal_axis == 1)
      {
        vl.x = vg.z;
        vl.y = vg.x;
        vl.z = vg.y;
      }
      else if (slab_normal_axis == 2)
      {
        vl.x = vg.x;
        vl.y = vg.y;
        vl.z = vg.z;
      }
      return vl;
    };

    inline Vector<double> give_slab_global_coordinates(const Vector<double> &vl)
    {
      Vector<double> vg;
      if (slab_normal_axis == 0)
      {
        vg.x = vl.y;
        vg.y = vl.z;
        vg.z = vl.x;
      }
      else if (slab_normal_axis == 1)
      {
        vg.x = vl.z;
        vg.y = vl.x;
        vg.z = vl.y;
      }
      else if (slab_normal_axis == 2)
      {
        vg.x = vl.x;
        vg.y = vl.y;
        vg.z = vl.z;
      }
      return vg;
    };

    bool initialized, calculated_once;
    void initialize();
    void make_slab_k_vectors();
    void make_slab_chi_vectors();

    void calculate_dipole_sum();

    double dipole_energy();
    double dipole_potential(const Vector<double> &);
    double dipole_potential(const int);
    Vector<double> dipole_field();

    double k_electrostatic;
    Vector<double> external_field;

    // simulation box lengths and its product in slab local coordinates.
    // the local coordinates depends on the choice of 'slab_normal_axis'.
    double lx, ly, lz, hz, lx_ly_inv;
    double slab_sum_e_coef;
    int kx_max, ky_max;
    std::vector<double> kx, ky, kp;
    std::vector<double> kx_coef, ky_coef, kp_coef;

    double dipole_coef;
    Vector<double> dipole_sum;          // dipole_sum: Sum_j(q_j vec(r_j))
    Vector<double> dipole_field_vector; // it's different from ewald_k's dipole.
                                        // it comes from change in summation order
                                        // from spherical to cylindrical.

    std::vector<Vector<double>> k_vector;
    std::vector<double> k_vector_sq;
    std::vector<double> field_k_coef; //, potential_k_coef;

    std::vector<std::vector<std::vector<std::vector<double>>>> chi_p;

    std::vector<std::vector<std::vector<double>>> chi_x;

    std::vector<std::vector<std::vector<double>>> chi_y;
  };

} // force_field

CAVIAR_NAMESPACE_CLOSE

#endif
