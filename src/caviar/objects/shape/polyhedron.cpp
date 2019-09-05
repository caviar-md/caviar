
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

#include "caviar/objects/shape/polyhedron.h"
#include "caviar/objects/shape/polyhedron/handler.h"
#include "caviar/utility/interpreter_io_headers.h"

namespace caviar {
namespace objects {
namespace shape {

Polyhedron::Polyhedron (CAVIAR *fptr) : Shape {fptr},
  polyhedron_handler {new shape::polyhedron::Handler{fptr}}
  {
  FC_OBJECT_INITIALIZE_INFO
}

Polyhedron::~Polyhedron () { 
  delete polyhedron_handler;  
}

bool Polyhedron::read (caviar::interpreter::Parser * parser) {
  FC_OBJECT_READ_INFO
  return polyhedron_handler -> read (parser);
}

bool Polyhedron::is_inside(const Vector<double> &v) {
  return polyhedron_handler -> is_inside (v); 
}

bool Polyhedron::is_inside(const Vector<double> &v, const double r) {
  return polyhedron_handler -> is_inside (v, r); 
}


bool Polyhedron::in_contact(const Vector<double> &v, const double r, Vector<double> & contact_vector) {
  return polyhedron_handler -> in_contact(v, r, contact_vector);
}

} //shape
} //objects
} // namespace caviar


