
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

#ifndef CAVIAR_OBJECTS_CONSTRAINT_CMMOTION_H
#define CAVIAR_OBJECTS_CONSTRAINT_CMMOTION_H

#include "caviar/objects/constraint.h"

CAVIAR_NAMESPACE_OPEN


namespace constraint {

/**
 * This class monitors and removes center of mass (CM) motions if requested.
 * This overall CM motion may comes from applying thermostats  directly to the atomic velocities 
 * instead of internal velocities. This phenomenon and internal (peculiar) velocities are
 * introduced in 'Adv. Polym. Sci. (2005) 173:105–149 DOI:10.1007/b99427'.
 * In the current version, we use this remedy. TODO. we may change this solution (as an option)
 * in the future releases.
 */
class Cm_motion : public Constraint {
 public:
  Cm_motion (class CAVIAR *);
   ~Cm_motion ( );
  bool read (class caviar::interpreter::Parser *);

  void apply_on_velocity (int64_t);


  void fix_velocity();
  void fix_angular_momentum();

  int velocity_steps;
  int angular_momentum_steps;  

};

} //constraint

CAVIAR_NAMESPACE_CLOSE

#endif
