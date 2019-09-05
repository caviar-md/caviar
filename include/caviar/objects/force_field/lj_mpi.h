
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_LJMPI_H
#define CAVIAR_OBJECTS_FORCEFIELD_LJMPI_H

#include "caviar/objects/force_field.h"

namespace caviar {
namespace objects {

namespace force_field {

/**
 * This class calculates LJ potential for the particles in the mpi mode
 * 
 */
class Lj_mpi : public Force_field {
public:
  Lj_mpi (class CAVIAR *);
  ~Lj_mpi () {};
  
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:
  std::vector<std::vector<Real_t>> epsilon,sigma;

};

} //force_field
} //objects
} // namespace caviar

#endif
