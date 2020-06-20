
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

namespace caviar {
namespace objects {
class Integrator;
class Atom_data;

/**
 * This class is the base class for all the constraints.
 * A constraint is called during a simulations. It can be implemented in
 * different parts of the virtual functions to be called at the desired moments
 * see the md_simulator objects to understant the time of call.
 */
class Constraint : public Pointers {
 public:
  Constraint (class CAVIAR *);
  virtual ~Constraint ( );
  virtual bool read (class caviar::interpreter::Parser *) = 0;

  virtual void step_part_I (int);
  virtual void step_part_II (int);
  virtual void step_part_III (int);

  /**
   *  Each integrator has a integer type. This is provided in cases the type
   *  matters for the constraint integration.
   */
  int integrator_type;
  std::shared_ptr<class Integrator > integrator;

  std::shared_ptr<class objects::Atom_data > atom_data;
 public:

  FC_BASE_OBJECT_COMMON_TOOLS
};

} //objects

} // namespace caviar

#endif
