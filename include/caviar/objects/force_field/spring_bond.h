
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_SPRINGBOND_H
#define CAVIAR_OBJECTS_FORCEFIELD_SPRINGBOND_H

#include "caviar/objects/force_field.h"

CAVIAR_NAMESPACE_OPEN

namespace force_field {

/**
 * This class does a spring force-field on the molecular bonds
 *  
 */
class Spring_bond : public Force_field {
public:
  Spring_bond (class CAVIAR *);
  ~Spring_bond () {};

  //double energy();

  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:

  std::vector<Real_t> elastic_coef, dissip_coef;
 
};

} //force_field

CAVIAR_NAMESPACE_CLOSE

#endif
