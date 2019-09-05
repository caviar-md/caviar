
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

#include "caviar/objects/shape/polyhedron/input.h"
#include "caviar/objects/shape/polyhedron/format_vtk_reader.h"
#include "caviar/objects/shape/polyhedron/format_stl_reader.h"

#include <string>
#include <fstream>

namespace caviar {
namespace objects {
namespace shape {
namespace polyhedron {

Input::Input (CAVIAR *fptr) : Pointers{fptr} {}

Input::~Input () { }

void Input::read_vtk (shape::polyhedron::Polyhedron & p_object, const std::string &file_name) {
  class Format_vtk_reader fvr (fptr);
  fvr.read_polyhedron (p_object, file_name);
}

void Input::read_stl (shape::polyhedron::Polyhedron & p_object, const std::string &file_name) {
  class Format_stl_reader fvr (fptr);
  fvr.read_polyhedron (p_object, file_name);
}

} //polyhedron
} //shape
} //objects
} // namespace caviar

