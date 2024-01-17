
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

#include "caviar/objects/shape/plane.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace shape
{

  Plane::Plane(CAVIAR *fptr) : Shape{fptr}, flatness_tol{0.001} {
                                                FC_OBJECT_INITIALIZE_INFO} Plane::~Plane() {}

  bool Plane::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR
      ASSIGN_REAL_3D_VECTOR(center, "Plane read: ", "")
      else ASSIGN_REAL_3D_VECTOR(normal, "Plane read: ", "") else ASSIGN_REAL(radius, "Plane read:", "") else ASSIGN_REAL(flatness_tol, "Plane read:", "") else error->all(FC_FILE_LINE_FUNC_PARSE, " Plane read: Unknown variable or command");
    }

    return in_file;
    ;
  }

  bool Plane::on_the_plane(const Vector<double> &v)
  {
    Vector<double> v_1 = v - center;
    if (v_1 * normal > flatness_tol)
      return false;
    return true;
  }

  bool Plane::is_inside(const Vector<double> &v)
  {
    if (!on_the_plane(v))
      return false;
    Vector<double> v_1 = v - center;
    if (v_1 * v_1 > radius * radius)
      return false;
    return true;
  }

  bool Plane::is_inside(const Vector<double> &v, const double r)
  {
    if (!on_the_plane(v))
      return false;
    Vector<double> v_1 = v - center;
    if (v_1 * v_1 > (radius - r) * (radius - r))
      return false;
    return true;
  }

  bool Plane::in_contact(const Vector<double> &v, const double r, Vector<double> &contact_vector)
  {
    std::string s = "incomplete function:";
    s += __FILE__ + std::to_string(__LINE__) + __func__;
    output->warning(s);
    std::cout << "  " << v << r << contact_vector << std::endl;
    return false;
  }

} // shape

CAVIAR_NAMESPACE_CLOSE
