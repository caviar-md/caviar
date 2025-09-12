
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATICEWALD1D_H
#define CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATICEWALD1D_H

#include "caviar/objects/force_field.h"

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  /**
   * This class has a one-dimentional-periodic electrostatic ewald force-field.
   *
   *
   */
  class Electrostatic_ewald1d : public Force_field
  {
  public:
    Electrostatic_ewald1d(class CAVIAR *);
    ~Electrostatic_ewald1d(){};
    double potential(const Vector<double> &);
    double potential(const int);

    double potential_r(const Vector<double> &);
    double potential_r(const int);

    double potential_k(const Vector<double> &);
    double potential_k(const int);

    Vector<double> field(const Vector<double> &);
    Vector<double> field(const int);

    Vector<double> field_r(const Vector<double> &);
    Vector<double> field_r(const int);

    Vector<double> field_k(const Vector<double> &);
    Vector<double> field_k(const int);

    double energy();

    bool read(class caviar::interpreter::Parser *);
    void verify_settings();
    void calculate_acceleration();

  public:
    std::vector<std::vector<Real_t>> lambda;
    bool lambda_is_set = false;
    double k_electrostatic;
    double sigma; //   smoothing-out parameter
    std::vector<Vector<double>> lattice_vec;
    int num_mirrors;
  };

} // force_field

CAVIAR_NAMESPACE_CLOSE

#endif
