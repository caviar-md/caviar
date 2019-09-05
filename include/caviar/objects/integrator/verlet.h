
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

#ifndef CAVIAR_OBJECTS_INTEGRATOR_VERLET_H
#define CAVIAR_OBJECTS_INTEGRATOR_VERLET_H

#include "caviar/objects/integrator.h"

namespace caviar {
namespace objects {
class Domain;
namespace integrator {

/**
 * This class has verlet method of integration.
 * Note that the constraint method based on velocity-rescaling does not affect
 * this type of integration.
 * The update formula are, 
 *
 * step_part_I ():
 *
 *  \f$  r (t + dt) = 2*r(t) − r (t − dt) + a(t) * dt^2  \f$
 *
 *  \f$  v(t) = (r (t + dt) - r (t - dt) )/ (2*dt)  \f$
 *
 * step_part_II (): 
 *
 *   NOTHING
 *
 * and is done in 'step_part_I ()' function.
 * TODO. The initialization formula 
 * \f$  r (dt) ~= r i (0) + v (0)dt + 1/2 * a(0) * dt^2  \f$
 * is not implemented yet. This should not make much a problem if one has equilibrium
 * steps prior to the main simulations.
 *
 * XXX. note that this is 'v(t)' not 'v(t+dt)'
 * This is implemented according to.
 * 'Computational Physics - Molecular Dynamics Simulations-
 * E. Carlon, M. Laleman and S. Nomidis – Academic year 2015/2016'
 * and also,
 * 'An overview of integration schemes for molecular dynamics simulations-
 * Ulf D. Schiller - 5th March 2008'.
 *
 */
class Verlet : public Integrator {
public:
  Verlet (class CAVIAR *);
   ~Verlet();
  bool read (class caviar::interpreter::Parser *);

  void step_part_I ();
  void step_part_II ();
  void step_part_III ();
  void verify_settings ();

  class objects::Domain *domain;

public:

};

} //integrator
} //objects
} // namespace caviar

#endif
