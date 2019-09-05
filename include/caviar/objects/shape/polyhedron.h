
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

#ifndef CAVIAR_OBJECTS_SHAPE_POLYHEDRON_H
#define CAVIAR_OBJECTS_SHAPE_POLYHEDRON_H

#include "caviar/objects/shape.h"

namespace caviar {
namespace objects {
namespace shape {
namespace polyhedron { class Handler; }

/**
 * This class has a polyhedron shape.
 * the shape must be imported in a geomety format such as vtk or stl
 * 
 */
class Polyhedron : public Shape {
  public:
  Polyhedron (class CAVIAR *);
  ~Polyhedron ();
  
  bool read(class caviar::interpreter::Parser *);
  
  bool is_inside (const Vector<double> &v);
  bool is_inside (const Vector<double> &, const double rad);   
  bool in_contact (const Vector<double> &, const double rad, Vector<double> & contact_vector);    

  public:
  void command_parameters (class caviar::interpreter::Parser *);    
  void command_generate ();

  class shape::polyhedron::Handler * polyhedron_handler;
  bool polyhedron_add;

};
} //shape
} //objects
} // namespace caviar

#endif
