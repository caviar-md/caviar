
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_LJCELLLIST_H
#define CAVIAR_OBJECTS_FORCEFIELD_LJCELLLIST_H

#include "caviar/objects/force_field.h"

namespace caviar {

namespace force_field {

/**
 * This class calculates LJ potential for the particles. It specially designed
 * for a Cell-list objects.
 */
class Lj_cell_list : public Force_field {
public:
  Lj_cell_list (class CAVIAR *);
  ~Lj_cell_list () {};
  
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:
  std::vector<std::vector<Real_t>> epsilon,sigma;  
};

} //force_field

} // namespace caviar

#endif
