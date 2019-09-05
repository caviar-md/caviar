
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

#ifndef CAVIAR_OBJECTS_INTEGRATOR_VELOCITYVERLETLANGEVIN_H
#define CAVIAR_OBJECTS_INTEGRATOR_VELOCITYVERLETLANGEVIN_H

#include "caviar/objects/integrator.h"
#include <random>

namespace caviar {
namespace objects {
namespace integrator {

/**
 * This class has classical velocity-verlet (or Leapfrog) method with Langevin thermostat implicitly.
 * It is introduced in
 * 'M. Kröger, 'Langevin dynamics modified Velocity Verlet algorithm'
 *
 *
 * The 'Classical Velocity Verlet algorithm' is called as another 
 * implementaion of Leapfrog at Rapaport 'Art of Molecular Dynamics'. 
 *
 *   \f$  a = (2.0 - friction * dt) / (2.0 + friction * dt)  \f$
 *
 *   \f$  b = std::sqrt(kb  * temperature * friction * 0.5 * dt)  \f$
 *
 *   \f$  c = 2.0 * dt / (2.0 + friction * dt)  \f$
 *
 *  step_part_I ():
 * @code
 *    eta [i] = random_vector
 *  
 *    vel [i] += 0.5 * acc [i] * dt + b * eta;
 *
 *    pos [i] += vel [i] * c;
 *  
 * @endcode
 *  step_part_II ():
 *
 * <code>
 *    vel [i] = a * vel [i] + b * eta[i] + 0.5 * acc [i] * dt; </code>
 * 
 */
class Velocity_verlet_langevin : public Integrator {
public:
  Velocity_verlet_langevin (class CAVIAR *);
   ~Velocity_verlet_langevin();    
  bool read (class caviar::interpreter::Parser *);
public:

  void step_part_I ();
  void step_part_II ();
  void step_part_III ();
  void verify_settings ();
  void initialize ();
  void print_langevin_parameters();

  std::mt19937 rnd_generator_x, rnd_generator_y, rnd_generator_z;
  std::vector<double> eta_x, eta_y, eta_z;
// the default stddev is 1  
  std::normal_distribution<double> rnd_ndist_x, rnd_ndist_y, rnd_ndist_z;    
  double temperature, friction, kb, kbt;
  double a, b, c; 
  bool initialized;
};

} //integrator
} //objects
} // namespace caviar

#endif
