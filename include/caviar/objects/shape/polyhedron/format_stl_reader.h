
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

#ifndef CAVIAR_OBJECTS_SHAPE_POLYHEDRON_FORMAT_STL_READER_h
#define CAVIAR_OBJECTS_SHAPE_POLYHEDRON_FORMAT_STL_READER_h

#include "caviar/utility/pointers.h"

namespace caviar {

namespace shape {
namespace polyhedron {
struct Polyhedron;
class Format_stl_reader : public Pointers {
public:

  Format_stl_reader (class CAVIAR *);
  ~Format_stl_reader();
  
  void read_polyhedron (shape::polyhedron::Polyhedron &, const std::string &);

// It checks if the two vertices are similar 
// Then makes a map of all vertices to the similar
// ones with the lower index, then clear
// void merge_vertices (int); 

};
} //polyhedron
} //shape

} // namespace caviar

#endif
