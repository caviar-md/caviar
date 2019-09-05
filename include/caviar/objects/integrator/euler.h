
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

#ifndef CAVIAR_OBJECTS_INTEGRATOR_EULER_H
#define CAVIAR_OBJECTS_INTEGRATOR_EULER_H

#include "caviar/objects/integrator.h"

namespace caviar {
namespace objects {
namespace integrator {

/**
 * This class has Euler method of integration. It is not a good method because
 * it is not time-reversible or phase-space preserving. We have implemented it 
 * for educational reasons.
 *
 * step_I ():
 *
 *  \f$ r (t+dt) = r(t) + v(t)*dt + 1/2 * a(t) * dt^2 \f$
 *
 *  \f$ v (t+dt) = v(t) + a(t)*dt \f$
 *
 * step_II ():
 *
 * NOTHING
 */
class Euler : public Integrator {
public:
  Euler (class CAVIAR *);
   ~Euler();
  bool read (class caviar::interpreter::Parser *);

  void step_part_I ();
  void step_part_II ();
  void step_part_III ();
  void verify_settings ();
public:

};

} //integrator
} //objects
} // namespace caviar

#endif
