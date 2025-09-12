
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

#pragma once

#include "caviar/objects/shape.hpp"

namespace caviar
{

  namespace shape
  {

    /**
     * This class has a plane shape.
     *
     *
     */
    class Plane : public Shape
    {
    public:
      Plane(class CAVIAR *);
      ~Plane();

      // there are different ways to define a circle: 3 points on it, centre and one point on it, centre and radius and normal,...
      bool read(class caviar::interpreter::Parser *);
      double radius;
      double flatness_tol;
      Vector3d<double> center;
      Vector3d<double> normal;
      bool on_the_plane(const Vector3d<double> &v);
      bool is_inside(const Vector3d<double> &v);
      bool is_inside(const Vector3d<double> &, const double rad);
      bool in_contact(const Vector3d<double> &, const double rad, Vector3d<double> &contact_vector);

      bool make_basis_vectors();
    };

  } // shape

}
