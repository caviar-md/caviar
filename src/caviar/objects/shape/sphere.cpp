
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

#include "caviar/objects/shape/sphere.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace shape
{

  Sphere::Sphere(CAVIAR *fptr) : Shape{fptr} {
                                     FC_OBJECT_INITIALIZE_INFO} Sphere::~Sphere() {}

  bool Sphere::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR
      ASSIGN_REAL_3D_VECTOR(center, "SPHERE read: ", "")
      else ASSIGN_REAL(radius, "SPHERE read:", "") else error->all(FC_FILE_LINE_FUNC_PARSE, " SPHERE read: Unknown variable or command");
    }

    return in_file;
    ;
  }

  bool Sphere::is_inside(const Vector<double> &v)
  {
    Vector<double> v_1 = v - center;
    if (v_1 * v_1 > radius * radius)
      return false;
    return true;
  }

  bool Sphere::is_inside(const Vector<double> &v, const double r)
  {
    Vector<double> v_1 = v - center;
    if (v_1 * v_1 > (radius - r) * (radius - r))
      return false;
    return true;
  }

  bool Sphere::in_contact(const Vector<double> &v, const double r, Vector<double> &contact_vector)
  {
    Vector<double> v_dif = v - center;
    const auto v_dif_sq = v_dif * v_dif;
    if (v_dif_sq < (radius - r) * (radius - r))
      return false;
    if (v_dif_sq > (radius + r) * (radius + r))
      return false;
    const auto v_dif_sq_norm = std::sqrt(v_dif_sq);
    const auto v_dif_norm = v_dif / v_dif_sq_norm;
    Vector<double> tmp = -v_dif_norm * (v_dif_sq_norm - radius);
    contact_vector += tmp;
    return true;
  }

} // shape

CAVIAR_NAMESPACE_CLOSE
