
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

#ifndef CAVIAR_OBJECTS_CONSTRAINT_H
#define CAVIAR_OBJECTS_CONSTRAINT_H

#include "caviar/utility/objects_common_headers.h"

CAVIAR_NAMESPACE_OPEN

class Integrator;
class Atom_data;

enum class Constraint_t
{
  Atom_molarity,
  Atoms_molarity,
  Berendsen,
  Cm_motion,
  M_shake,
  Nose_hoover,
  Nve,
  Rattle,
  Shake,
  Unknown
};

/**
 * This class is the base class for all the constraints.
 * A constraint is called during a simulations. It can be implemented in
 * different parts of the virtual functions to be called at the desired moments
 * see the md_simulator objects to understant the time of call.
 */
class Constraint : public Pointers
{
public:
  Constraint(class CAVIAR *);
  virtual ~Constraint();
  virtual bool read(class caviar::interpreter::Parser *) = 0;

/**
*   it will be applied at the end of each time step
*/
  virtual void apply(int64_t);


/**
 * it should be used after calculating new position to fix it before using it for acceleration calculation.
 * for example, shake, m-shake, rattle use it
*/
  virtual void fix_position(int64_t); 

/**
 * it should be used after calculating velocity to fix it. 
*/
  virtual void fix_velocity(int64_t, bool &recalculate_temperature); 

/**
 * it should be used after calculating acceleration to fix it.
*/  
virtual void fix_acceleration(int64_t); 
  

/**
 * shake, m-shake, rattle use it
*/
  virtual void apply_shake(int64_t); 

/**
 * NVE, NVT fix. It applies on velocity, but it must be called only on specific algorithm step
*/
  virtual void apply_thermostat(int64_t, bool &recalculate_temperature); 

/**
 * NPT fix. It applies on position, but it must be called only on specific algorithm step.
 * if (fix_position_needed) returns true, it means that the barostat scaling is applied
 * on position and one needs to call shake algorithm again
*/
  virtual void apply_barostat(int64_t, bool &fix_position_needed); 


  class Atom_data *atom_data = nullptr;

  Constraint_t constraint_type;

  /**
   * MPI rank of the classs
  */
  int my_mpi_rank = -1;
public:
  FC_BASE_OBJECT_COMMON_TOOLS
};

CAVIAR_NAMESPACE_CLOSE

#endif
