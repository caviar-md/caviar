
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
     * This class has a polygon shape.
     *
     *
     */
    class Polygon : public Shape
    {
    public:
      Polygon(class CAVIAR *);
      ~Polygon();

      bool read(class caviar::interpreter::Parser *);
      // bool read (caviar::interpreter::Parser *, class Object_container *);

      bool is_inside(const Vector3d<double> &v);
      bool is_inside(const Vector3d<double> &, const double rad);
      bool in_contact(const Vector3d<double> &, const double rad, Vector3d<double> &contact_vector);

      std::vector<Vector3d<double>> p_3D;            // vertex points of the shape in 3D space
      Vector3d<double> n_vector, u_vector, v_vector; //  basis vectors of the plane: normal and (u,v)
      std::vector<double> u_list, v_list;          // array containing u and v coordinates of the vertices in the plane
      double mat_inv[3][3];                        // inverse of the transformation matrix;

      double flatness_tol;

      void make_basis_vectors();
      void make_uv_vectors();
      void make_transform_matrix();
      void uv_to_xyz(double u_i, double v_i, Vector3d<double> &p_i);
      void xyz_to_uv(const Vector3d<double> &p_i, double &u_i, double &v_i);
      bool is_inside(double u_i, double v_i); // using Jordan curve theorem,
                                              // checks whether the point is inside the curve or not;
                                              // it has a bug when (v_list[i] == v_list[j]);
    };

  } // shape

}
