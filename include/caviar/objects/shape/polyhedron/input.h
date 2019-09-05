
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

#ifndef CAVIAR_OBJECTS_SHAPE_POLYHEDRON_INPUT_H
#define CAVIAR_OBJECTS_SHAPE_POLYHEDRON_INPUT_H

#include "caviar/utility/pointers.h"

namespace caviar {
namespace objects {
namespace shape {
namespace polyhedron {
struct Polyhedron;
class Input : public Pointers {
public:
  Input (class CAVIAR *);
  ~Input ();
  
  void read_vtk (shape::polyhedron::Polyhedron&, const std::string &); 
  void read_stl (shape::polyhedron::Polyhedron&, const std::string &); 


};
} //polyhedron
} //shape
} //objects
} // namespace caviar

#endif
