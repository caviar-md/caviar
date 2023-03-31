
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_FIXBOND_H
#define CAVIAR_OBJECTS_FORCEFIELD_FIXBOND_H

#include "caviar/objects/force_field.h"

CAVIAR_NAMESPACE_OPEN

namespace force_field {

/**
 * This class adds a atomic bond between atoms/molecules
 * if the given conditions are satisfied.
 */
class Fix_bond : public Force_field {
public:
  Fix_bond (class CAVIAR *);
  ~Fix_bond () {};


  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
  void create_atomic_bond ();
public:

  int type_i;
  int type_j;
  int btype;
  double Rmin;
  double blength;
  
  /**
   * Maximum number of bonds allowed for a atom type.
   * This function uses 'Atom_data.owned.atomic_bond_count' vector.
   */
  std::vector<int> bond_limit;
 
};

} //force_field

CAVIAR_NAMESPACE_CLOSE

#endif
