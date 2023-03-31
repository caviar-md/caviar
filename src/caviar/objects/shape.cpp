
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

#include "caviar/objects/shape.h"

CAVIAR_NAMESPACE_OPEN



Shape::Shape (CAVIAR *fptr) : Pointers{fptr} {
  FC_OBJECT_INITIALIZE
}

Shape::~Shape () {}

void Shape::verify_settings () {
  
}

bool Shape::is_outside (const Vector<double> &v) {
  return !is_inside(v);
}
  
bool Shape::is_outside (const Vector<double> &v, const double r) {
  return !is_inside(v, r);
}


CAVIAR_NAMESPACE_CLOSE

