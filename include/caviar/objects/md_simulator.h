
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

namespace caviar {

namespace objects {
class Atom_data;
class Integrator;
class Neighborlist;
class Force_field;
class Constraint;
class Writer;

/**
 * This class is the base class for all the md_simulators.
 * It relates the required elements of a MD simulation.
 * 
 */
class Md_simulator : public Pointers {
 public:
  Md_simulator (class CAVIAR *);
  virtual ~Md_simulator ( );
  virtual bool read (class caviar::interpreter::Parser *) = 0;
  virtual bool read_base_class_commands(class caviar::interpreter::Parser *);
  virtual bool run () = 0;
  virtual void step(int i);
  virtual void initialize();
  virtual void setup ();
  virtual void cleanup ();
  virtual bool boundary_condition();

  class objects::Atom_data *atom_data;
  class objects::Integrator *integrator;
  std::vector<class objects::Neighborlist *> neighborlist;
  std::vector<objects::Force_field *> force_field; 
  std::vector<objects::Constraint *> constraint; 
  std::vector<objects::Writer *> writer; 

  double time, dt;
  clock_t t_start, t_end;
  double initial_time, final_time;
  bool use_time, use_step;
  int initial_step, final_step;
  bool initialized;

  FC_BASE_OBJECT_COMMON_TOOLS
};

} //objects

} // namespace caviar

#endif
