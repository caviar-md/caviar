
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

#ifndef CAVIAR_OBJECTS_SHAPE_POLYHEDRON_POINT_INSIDE_H
#define CAVIAR_OBJECTS_SHAPE_POLYHEDRON_POINT_INSIDE_H

#include "caviar/utility/pointers.h"
#include "caviar/utility/vector.h"

CAVIAR_NAMESPACE_OPEN

namespace shape {
namespace polyhedron {
struct Polyhedron;

/**
 * A class to understant if a point is inside the polyhedron or if a particle is
 * in contact with it.
 */
class Point_Inside : public Pointers {
public:
  Point_Inside (class CAVIAR *);
  ~Point_Inside ();

  /**
   *
   * usage in forcefields
   * 'in_contact_grid()' and 'in_contact_all()' are both do the same thing, and
   * should give the same results.
   * except one calculates the distance to all the polygons and the other do it
   * only to the polygons that are in the grid.
   * the grid one may be faster (depending on the number of the polygons) but prone to
   * bugs if the parameters are not correct.
   */
  bool in_contact_grid (shape::polyhedron::Polyhedron &, const Vector<Real_t> &v1, const Real_t radius, Vector<double> &contact_vector);
  bool in_contact_all  (shape::polyhedron::Polyhedron &, const Vector<Real_t> &v1, const Real_t radius, Vector<double> &contact_vector);


  /**  
   * usage in particle initial distributions (random, grid or ...)
   */
  bool is_inside(shape::polyhedron::Polyhedron &, const Vector<double> &v0);

  /** 
   * usage in particle initial distributions (random, grid or ...)
   */
  bool is_inside_grid(shape::polyhedron::Polyhedron &, const Vector<double> &v, const double r) ;
  bool is_inside_all(shape::polyhedron::Polyhedron &, const Vector<double> &v, const double r) ;

  /**
   * usage in particle initial distributions (random, grid or ...)
   * when we just to know a sphere with radius r is contacted the surface and so
   * it is not permitted to be created.
   */
  bool in_contact_grid (shape::polyhedron::Polyhedron &, const Vector<Real_t> &v1, const Real_t radius);
  bool in_contact_all  (shape::polyhedron::Polyhedron &, const Vector<Real_t> &v1, const Real_t radius);

  /**
   * usage in particle initial distributions (random, grid or ...)
   * as the name tells, by using a famous method that I don't recall its name now,
   * checks if a point is inside polyhedron by using an imaginary beam and counting
   * how many times it crosses polygons of the polyhedron.
   */
  bool ray_tells_point_is_inside (shape::polyhedron::Polyhedron &, const Vector<Real_t> &v1, const int ray_axis);
  int point_is_inside_method;

};
} //polyhedron
} //shape


CAVIAR_NAMESPACE_CLOSE

#endif
