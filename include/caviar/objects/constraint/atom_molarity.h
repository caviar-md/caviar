
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

#ifndef CAVIAR_OBJECTS_CONSTRAINT_ATOMMOLARITY_H
#define CAVIAR_OBJECTS_CONSTRAINT_ATOMMOLARITY_H

#include "caviar/objects/constraint.h"

namespace caviar {

namespace unique {
class Atom; class Molecule;
}
namespace constraint {

/**
 * This class fixes number of a single type of atom in an area
 * 
 * 
 */
class Atom_molarity : public Constraint {
 public:
  Atom_molarity (class CAVIAR *);
   ~Atom_molarity ( );
  bool read (class caviar::interpreter::Parser *);

  void apply (int64_t);

  void verify_settings();

  bool minimum_set, maximum_set;
  int maximum_limit;
  int atom_type;
  int minimum_limit;
  int creation_try;
  int steps, check_steps;
  Vector<double> calculation_box_low, calculation_box_high;
  Vector<double> creation_box_low,    creation_box_high;
  unique::Molecule *creation_molecule;
  unique::Atom *creation_atom;
  bool settings_verified;  

};

} //constraint

} // namespace caviar

#endif
