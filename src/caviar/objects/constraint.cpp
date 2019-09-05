
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

#include "caviar/objects/constraint.h"

namespace caviar {

namespace objects {

Constraint::Constraint (CAVIAR *fptr) : Pointers{fptr}, integrator_type{-1},
    integrator{nullptr}, atom_data{nullptr} {
  FC_OBJECT_INITIALIZE
}

Constraint::~Constraint () {}

void Constraint::verify_settings () {
  
}

void Constraint::step_part_I (int) { }
void Constraint::step_part_II (int) { }
void Constraint::step_part_III (int) { }

} //objects


} // namespace caviar

