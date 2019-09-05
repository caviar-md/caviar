
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATICSHORTRANGE_H
#define CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATICSHORTRANGE_H

#include "caviar/objects/force_field.h"

namespace caviar {
namespace objects {
namespace force_field {

/**
 * This class has a short-range type of electrostatic force-field.
 * introduced in the paper (Steinbach 1994)
 * journal of computational chemistry Vol 15 No 7 pages 667-683 (1994)
 *  V(r) = ( q_i q_j / 4 Pi epsilon ) * (1/r + c r ^beta + d) 
 */
class Electrostatic_short_range : public Force_field {
public:
  Electrostatic_short_range (class CAVIAR *);
  ~Electrostatic_short_range () {};
  double potential (const Vector<double> &);
  double potential (const int);

  Vector<double> field (const Vector<double> &);
  Vector<double> field (const int);

  double energy();

  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();

  void initialize();
  double k_electrostatic;
  Vector<double> external_field;

  // short range force shift
  // V(r) = ( q_i q_j / 4 Pi epsilon ) * (1/r + c r ^beta + d) 
  // (Steinbach 1994)
  // journal of computational chemistry Vol 15 No 7 pages 667-683 (1994)
  double C, D, beta; 
  bool initialized;
  double cutoff_sq;
 
};

} //force_field
} //objects
} // namespace caviar

#endif
