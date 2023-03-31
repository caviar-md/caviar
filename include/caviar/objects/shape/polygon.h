
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

#ifndef CAVIAR_OBJECTS_SHAPE_POLYGON_H
#define CAVIAR_OBJECTS_SHAPE_POLYGON_H

#include "caviar/objects/shape.h"

CAVIAR_NAMESPACE_OPEN

namespace shape {

/**
 * This class has a polygon shape.
 * 
 * 
 */
class Polygon : public Shape {
public:
  Polygon (class CAVIAR *) ;
  ~Polygon ();  
  
  bool read(class caviar::interpreter::Parser *);
  //bool read (caviar::interpreter::Parser *, class Object_container *);  

  bool is_inside (const Vector<double> &v);
  bool is_inside (const Vector<double> &, const double rad);  
  bool in_contact (const Vector<double> &, const double rad, Vector<double> & contact_vector);
  
  std::vector<Vector<Real_t>> p_3D; // vertex points of the shape in 3D space
  Vector<Real_t> n_vector, u_vector, v_vector; //  basis vectors of the plane: normal and (u,v)
  std::vector<Real_t> u_list, v_list; // array containing u and v coordinates of the vertices in the plane
  Real_t mat_inv[3][3]; // inverse of the transformation matrix;

  double flatness_tol;
  

  void make_basis_vectors() ;
  void make_uv_vectors() ; 
  void make_transform_matrix();
  void uv_to_xyz (Real_t u_i, Real_t v_i, Vector<Real_t> &p_i);
  void xyz_to_uv (const Vector<Real_t> &p_i, Real_t &u_i, Real_t &v_i) ;
  bool is_inside (double u_i, double v_i); // using Jordan curve theorem,
      // checks whether the point is inside the curve or not;
      // it has a bug when (v_list[i] == v_list[j]);

};

} //shape

CAVIAR_NAMESPACE_CLOSE

#endif
