
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

#ifndef CAVIAR_OBJECTS_SHAPE_POLYHEDRON_POSTPROCESS_H
#define CAVIAR_OBJECTS_SHAPE_POLYHEDRON_POSTPROCESS_H

#include "caviar/utility/pointers.h"

namespace caviar {

namespace shape {
namespace polyhedron {
struct Polyhedron;
class Postprocess : public Pointers{
public:
  Postprocess (class CAVIAR *);
  ~Postprocess ();
  
  /**
   * contains the faces neccesary to check
   * make_grid has to be called after lowest_highest_coord()
   */
  void make_grid (shape::polyhedron::Polyhedron &); 

  /**
   * calculates gxlo, gxhi, gylo...
   */
  void lowest_highest_coord (shape::polyhedron::Polyhedron &); 


};
} //polyhedron
} //shape


} // namespace caviar

#endif
