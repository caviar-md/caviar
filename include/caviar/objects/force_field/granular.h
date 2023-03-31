
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_GRANULAR_H
#define CAVIAR_OBJECTS_FORCEFIELD_GRANULAR_H

#include "caviar/objects/force_field.h"


namespace caviar {

namespace force_field {

/**
 * This class makes simple Hookean force-field for the particles;
 * force_type = 0; linear dashpot
 * force_type = 1; visco-elastic
 */
class Granular : public Force_field {
public:
  Granular (class CAVIAR *);
  ~Granular () {};
  
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:
  std::vector<Real_t> elastic_coef, dissip_coef, radius;
  Vector<Real_t> gravity;
  int force_type; 
};

} //force_field

} // namespace caviar

#endif
