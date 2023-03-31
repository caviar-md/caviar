
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_DPD_H
#define CAVIAR_OBJECTS_FORCEFIELD_DPD_H

#include "caviar/objects/force_field.h"

#include <random>

/**
 * This class has dissipative particle dynamics force-field.
 * 
 * 
 */
CAVIAR_NAMESPACE_OPEN

namespace force_field {

class Dpd : public Force_field { // there's a numeric error due to using ghost atoms and different random number generated for their owned counterpart atom. // this problem is solved at Dpd_acc
public:
  Dpd (class CAVIAR *);
  ~Dpd () {};
  
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:
  std::vector<std::vector<Real_t>> conserv_coef,dissip_coef;
  Real_t temperature, kBoltzman;
  int rnd_seed;
  std::mt19937 rnd_generator;
  std::normal_distribution<double> rnd_ndist; // stddev() == 1
  double dt;
};

} //force_field

CAVIAR_NAMESPACE_CLOSE

#endif
