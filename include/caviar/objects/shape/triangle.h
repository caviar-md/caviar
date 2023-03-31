
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

#ifndef CAVIAR_OBJECTS_SHAPE_TRIANGLE_H
#define CAVIAR_OBJECTS_SHAPE_TRIANGLE_H

#include "caviar/objects/shape.h"

namespace caviar {

namespace shape {

/**
 * This class has a triangle shape. It could be used in order to make more 
 * complex shapes
 * 
 */
class Triangle : public Shape {
public:
  Triangle (class CAVIAR *) ;
  ~Triangle ();

// there are different ways to define a circle: 3 points on it, centre and one point on it, centre and radius and normal,...
  bool read(class caviar::interpreter::Parser *);
  double radius;
  double flatness_tol;
  Vector<double> center;
  Vector<double> normal;
  bool on_the_plane (const Vector<double> &v);
  bool is_inside (const Vector<double> &v);
  bool is_inside (const Vector<double> &, const double rad);  
  bool in_contact (const Vector<double> &, const double rad, Vector<double> & contact_vector);  
  
  bool make_basis_vectors();


};

} //shape

} // namespace caviar

#endif
