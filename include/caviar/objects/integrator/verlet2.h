
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

#ifndef CAVIAR_OBJECTS_INTEGRATOR_VERLET2_H
#define CAVIAR_OBJECTS_INTEGRATOR_VERLET2_H

#include "caviar/objects/integrator.h"

namespace caviar {
namespace objects {
namespace integrator {

/**
 * This class has an implementation of verlet method integration.
 * It is experimental and may be similar to other classes.
 * XXX.
 */
class Verlet2 : public Integrator {
public:
  Verlet2 (class CAVIAR *);
   ~Verlet2();
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
