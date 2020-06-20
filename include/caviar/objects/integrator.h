
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

#ifndef CAVIAR_OBJECTS_INTEGRATOR_H
#define CAVIAR_OBJECTS_INTEGRATOR_H

#include "caviar/utility/objects_common_headers.h"

namespace caviar {

namespace objects {
class Atom_data;

/**
 * This class is the base class for all the integrators for a MD.
 * 
 * 
 */
class Integrator : public Pointers {
public:
  Integrator (class CAVIAR *);
  virtual ~Integrator  ();
  virtual bool read (class caviar::interpreter::Parser *) = 0;

public:
  virtual void step_part_I () = 0;
  virtual void step_part_II () = 0;
  virtual void step_part_III () = 0;
  std::shared_ptr<class objects::Atom_data> atom_data;

  Real_t dt;

  // to be used in other objects;
  int integrator_type;

  FC_BASE_OBJECT_COMMON_TOOLS  
};

} //objects

} // namespace caviar

#endif
