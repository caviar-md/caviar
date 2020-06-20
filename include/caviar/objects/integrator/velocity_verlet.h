
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

#ifndef CAVIAR_OBJECTS_INTEGRATOR_VELOCITYVERLET_H
#define CAVIAR_OBJECTS_INTEGRATOR_VELOCITYVERLET_H

#include "caviar/objects/integrator.h"
//#include "caviar/utility/python_utils_dec.h"

namespace caviar {
namespace objects {
namespace integrator {

/**
 * This class has a velocity-verlet integration method.
 * It is implemented according to.
 * 'Computational Physics - Molecular Dynamics Simulations-
 * E. Carlon, M. Laleman and S. Nomidis – Academic year 2015/2016'
 * and also,
 * 'An overview of integration schemes for molecular dynamics simulations-
 * Ulf D. Schiller - 5th March 2008'.
 * We have implemented this algorithm in two types. The default is type 1,
 * which is the same as 'leap_frog.h' implementation. 
 *
 * TYPE 1 :  
 *
 *  step_part_I ():
 *
 *  \f$ r(t+dt) = r(t) + v(t)*dt + 1/2 * a(t) * dt^2 \f$
 *
 *  \f$ v (t + dt/2) = v (t) + (dt/2) a (t) \f$
 *
 *
 *
 *  step_part_II ():
 *
 *  \f$  v (t + dt) = v (t + dt/2) + (dt/2) a (t + dt) \f$
 *
 * TYPE 2 :  
 * 
 *  step_part_I ():
 *
 *  \f$  r(t+dt) = r(t) + v(t)*dt + 1/2 * a(t) * dt^2 \f$
 *
 *  step_part_II ():
 *
 *  \f$  v(t+dt) = v(t) + ( a(t+dt) + a(t) ) * dt / 2  \f$
 *
 * 
 */
class Velocity_verlet : public Integrator {
public:
  Velocity_verlet (class CAVIAR *);
   ~Velocity_verlet();
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
public:

  void step_part_I ();
  void step_part_II ();
  void step_part_III ();
  int type;
   
};

//void export_py_Velocity_verlet ();

} //integrator
} //objects
} // namespace caviar

#endif
