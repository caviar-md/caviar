
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

#ifndef CAVIAR_OBJECTS_INTEGRATOR_LEAPFROG2_H
#define CAVIAR_OBJECTS_INTEGRATOR_LEAPFROG2_H

#include "caviar/objects/integrator.h"

namespace caviar {
namespace objects {
namespace integrator {

/**
 * This class has another method of Leap-Frog implementation
 * The method is intoduced in  Rapaport 'Art of Molecular Dynamics'.
 * . XXX. It is in experimental stage.
 *
 *  step_part_I ():
 *
 *  \f$ v (t + dt/2) = v  (t − dt/2) + dt * a  (t)  \f$
 *
 *  \f$ r (t + dt) = r (t) + dt * v (t + dt/2) \f$
 *
 *  step_part_II ():
 *
 *  nothing
 * 
 */
class Leap_frog2 : public Integrator {
public:
  Leap_frog2 (class CAVIAR *);
   ~Leap_frog2();
  bool read (class caviar::interpreter::Parser *);

  void step_part_I ();
  void step_part_II ();
  void step_part_III ();
  void verify_settings ();
public:
  std::vector <Vector<double>> vel_h2;
};

} //integrator
} //objects
} // namespace caviar

#endif
