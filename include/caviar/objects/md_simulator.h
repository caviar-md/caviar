
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

#ifndef CAVIAR_OBJECTS_MDSIMULATOR_H
#define CAVIAR_OBJECTS_MDSIMULATOR_H

#include "caviar/utility/objects_common_headers.h"

#include <random> // used for velocity_verlet_langevin

namespace caviar {

namespace objects {
namespace unique
{
class Time_function;
class Time_function_3d;
}
class Atom_data;
class Integrator;
class Neighborlist;
class Force_field;
class Constraint;
class Writer;

enum class Integrator_t {
  Verlet,
  Velocity_verlet,
  Leap_frog,
  Velocity_verlet_langevin,
  Unknown
};

/**
 * This class is the base class for all the md_simulators.
 * It relates the required elements of a MD simulation.
 * 
 * ------------------------------------------------------------------
 * 
 * The algorithm with Velocity_verlet_langevin integration comes from
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
 *  Part I:
 * @code
 *    eta [i] = random_vector
 *  
 *    vel [i] += 0.5 * acc [i] * dt + b * eta;
 *
 *    pos [i] += vel [i] * c;
 *  
 * @endcode
 *  Part II:
 *
 * <code>
 *    vel [i] = a * vel [i] + b * eta[i] + 0.5 * acc [i] * dt; </code>
 *  
 * ----------------------------------
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
 * --------------------------------------
 * 
 * This class do Leap-Frog integration
 * 
 * The method is intoduced and implemented in  Rapaport 'Art of Molecular Dynamics'.
 *
 * It is the same as 'Classical Velocity Verlet algorithm' according to
 * M. Kröger, 'Langevin dynamics modified Velocity Verlet algorithm'
 * in which cited, M. Kröger,
 * Models for polymeric and anisotropic liquids
 * (Springer, Berlin, 2005).
 *
 *  step_part_I ():
 *
 *  \f$ v (t + dt/2) = v (t) + (dt/2) a (t) \f$
 *
 *  \f$ r (t + dt) = r (t) + dt * v (t + dt/2) \f$
 *
 *  step_part_II ():
 *
 *  \f$ v (t + dt) = v (t + dt/2) + (dt/2) a (t + dt) \f$
 * 
 * --------------------------------------
 * 
 * * This class has a velocity-verlet integration method.
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
 * -------------------------------------------------------------
 * 
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
class Md_simulator : public Pointers {
 public:
  Md_simulator (class CAVIAR *);
  virtual ~Md_simulator ( );
  virtual bool read (class caviar::interpreter::Parser *);
  virtual bool read_base_class_commands(class caviar::interpreter::Parser *);
  virtual bool run ();
  virtual void step();
  virtual void step(int64_t);
  virtual void initialize();
  virtual void setup ();
  virtual void cleanup ();
  virtual bool boundary_condition();
  //virtual void verify_settings ();
  virtual void re_calculate_acc();
  
  class objects::Atom_data *atom_data;
  //class objects::Integrator *integrator;
  std::vector<class objects::Neighborlist *> neighborlist;
  std::vector<objects::Force_field *> force_field; 
  std::vector<objects::Constraint *> constraint; 
  std::vector<objects::Writer *> writer; 
  std::vector<objects::unique::Time_function *> time_function; 
  std::vector<objects::unique::Time_function_3d *> time_function_3d; 
  
  double time, dt;
  clock_t t_start, t_end;
  double initial_time, final_time;
  bool use_time, use_step;
  int64_t current_step, initial_step, final_step;
  bool initialized;
  
  
  
  Integrator_t integrator_type;

  void integrate_leap_frog();
  void integrate_velocity_verlet();
  void integrate_velocity_verlet_langevin();

  struct {
    std::mt19937 rnd_generator_x, rnd_generator_y, rnd_generator_z;
    std::vector<double> eta_x, eta_y, eta_z;
  // the default stddev is 1  
    std::normal_distribution<double> rnd_ndist_x, rnd_ndist_y, rnd_ndist_z;    
    double temperature, friction, kb, kbt;
    double a, b, c; 
    
  } langevin_param;

  FC_BASE_OBJECT_COMMON_TOOLS
};

} //objects

} // namespace caviar

#endif
