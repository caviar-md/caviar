
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

#include "caviar/objects/force_field.h"
#include "caviar/interpreter/error.h"

namespace caviar {



Force_field::Force_field (CAVIAR *fptr) : Pointers{fptr}, 
  atom_data{nullptr}, domain{nullptr}, neighborlist{nullptr} {
  FC_OBJECT_INITIALIZE
}

void Force_field::verify_settings () {
  
}

double Force_field::energy() {
  error->all(FC_FILE_LINE_FUNC, "The energy calculation of this force_field is not implemented");
  return 0.0;
}

double Force_field::potential (const Vector<double> &) {
  error->all(FC_FILE_LINE_FUNC, "The potential calculation of this force_field is not implemented");
  return 0.0;
}

double Force_field::potential (const int) {
  error->all(FC_FILE_LINE_FUNC, "The potential calculation of this force_field is not implemented");
  return 0.0;
}


Vector<double> Force_field::field (const Vector<double> &) {
  error->all(FC_FILE_LINE_FUNC, "The field calculation of this force_field is not implemented");
  return Vector<double> {0,0,0};
}

Vector<double> Force_field::field (const int) {
  error->all(FC_FILE_LINE_FUNC, "The field calculation of this force_field is not implemented");
  return Vector<double> {0,0,0};
}




} // namespace caviar

