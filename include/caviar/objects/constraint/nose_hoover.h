
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

#ifndef CAVIAR_OBJECTS_CONSTRAINT_NOSEHOOVER_H
#define CAVIAR_OBJECTS_CONSTRAINT_NOSEHOOVER_H

#include "caviar/objects/constraint.h"

namespace caviar {

namespace constraint {

/**
 * This class has Nose-Hoover thermostat.
 * It is implemented according to
 * 'Berendsen and Nose-Hoover thermostats Victor Ruhle August 8, 2007'
 * 
 */
class Nose_hoover : public Constraint {
 public:
  Nose_hoover (class CAVIAR *);
  ~Nose_hoover ( );
  bool read (class caviar::interpreter::Parser *);

  void apply_on_acceleration (int64_t);

  void verify_settings();

  /**
   * This fictious mass determines the coupling between heat bath and the system.
   * Used in type 1 .
   */
  double mass;

  /**
   * effective relaxation time. Used in type 2 .
   */
  double tau;

  /**
   * additional degree of freedom related to the heat bath
   */
  double zeta, zeta_dot;

  /**
   * timestep value
   */
  double dt;

  bool settings_verified;  

  /**
   * the type of thermostat implementation
   */
  int type;

  /**
   * Boltzman constant. Used in type 1 .
   */
  double kb;

  double temperature;

 public:

};

} //constraint

} // namespace caviar

#endif
