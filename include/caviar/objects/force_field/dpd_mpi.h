
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_DPDMPI_H
#define CAVIAR_OBJECTS_FORCEFIELD_DPDMPI_H

#include "caviar/objects/force_field.h"

#include <random>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  /**
   * This class dpd force-field with mpi implementation
   *
   *
   */
  class Dpd_mpi : public Force_field
  { // this is maybe slower than dpd but more accurate, in the scense that it will send the force caclulated to the owned counterpart of ghost atoms.
  public:
    Dpd_mpi(class CAVIAR *);
    ~Dpd_mpi(){};

    bool read(class caviar::interpreter::Parser *);
    void verify_settings();
    void calculate_acceleration();

  public:
    std::vector<std::vector<Real_t>> conserv_coef, dissip_coef;
    Real_t temperature, kBoltzman;
    int rnd_seed;
    std::mt19937 rnd_generator;
    std::normal_distribution<double> rnd_ndist; // stddev() == 1
    double dt;
  };

} // force_field

CAVIAR_NAMESPACE_CLOSE

#endif
