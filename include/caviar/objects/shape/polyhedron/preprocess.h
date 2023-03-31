
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

#ifndef CAVIAR_OBJECTS_SHAPE_POLYHEDRON_PREPROCESS_H
#define CAVIAR_OBJECTS_SHAPE_POLYHEDRON_PREPROCESS_H

#include "caviar/utility/pointers.h"

CAVIAR_NAMESPACE_OPEN

namespace shape {
namespace polyhedron {

struct Polyhedron;
class Preprocess : public Pointers {
public:
  Preprocess (class CAVIAR *);
  ~Preprocess ();

  /**
   * checks neighborlist faces and sorts the vertices so that their normal vectors would be alighned when created.
   */
  void pre_correct_normals (shape::polyhedron::Polyhedron&); 

  void merge_vertices (shape::polyhedron::Polyhedron&);
  
};
} //polyhedron
} //shape

CAVIAR_NAMESPACE_CLOSE

#endif
