
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

#ifndef CAVIAR_OBJECTS_CONSTRAINT_BERENDSEN_H
#define CAVIAR_OBJECTS_CONSTRAINT_BERENDSEN_H

#include "caviar/objects/constraint.h"

namespace caviar {
namespace objects {
namespace constraint {

/**
 * This class has Brendsen thermostat implemented. 
 * It is implemented according to
 * 'Berendsen and Nose-Hoover thermostats Victor Ruhle August 8, 2007'
 * 
 */
class Berendsen : public Constraint {
 public:
  Berendsen (class CAVIAR *);
  ~Berendsen ( );
  bool read (class caviar::interpreter::Parser *);

  void step_part_I (int);

  void verify_settings();

  double coupling, dt, temperature;

 public:

};

} //constraint
} //objects
} // namespace caviar

#endif
