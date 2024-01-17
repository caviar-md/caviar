
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

#include "caviar/objects/shape/polygon.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace shape
{

  Polygon::Polygon(CAVIAR *fptr) : Shape{fptr},
                                   flatness_tol{0.001} {
                                       FC_OBJECT_INITIALIZE_INFO}

                                   Polygon::~Polygon()
  {
  }

  bool Polygon::is_inside(const Vector<double> &v)
  {
    // error->all (FC_FILE_LINE_FUNC_PARSE, "not_implemented");
    is_inside(v.x, v.y);
    output->warning("Polygon::is_inside: not implemented yet.");
    return true;
  }

  bool Polygon::is_inside(const Vector<double> &v, const double r)
  {
    // error->all (FC_FILE_LINE_FUNC_PARSE, "not_implemented");
    is_inside(v.x + r, v.y + r);
    output->warning("Polygon::is_inside: not implemented yet.");
    return true;
  }

  bool Polygon::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;
    std::vector<Vector<double>> tmp;
    int no_points = 0;
    while (true)
    {

      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR
      ++no_points;
      tmp.resize(no_points);
      // ASSIGN_REAL_3D_VECTOR(tmp[no_points-1], "Closed lines read: ","")
      GET_OR_CHOOSE_A_REAL_3D_VECTOR(tmp[no_points - 1], "Closed lines read: ", "")
      // else error->all (FC_FILE_LINE_FUNC_PARSE, " closed lines read: Unknown variable or command");
    }

    for (auto i : tmp)
      p_3D.push_back(i);

    return in_file;
  }

  void Polygon::make_basis_vectors()
  {
    if (p_3D.size() > 1)
    {
      Vector<Real_t> vec1 = p_3D[1] - p_3D[0];
      Vector<Real_t> vec2 = p_3D[2] - p_3D[0];
      std::vector<std::vector<Real_t>> mat_inv;
      n_vector = cross_product(vec1, vec2);
      u_vector = vec1;
      v_vector = cross_product(n_vector, u_vector);
      normalize(u_vector);
      normalize(v_vector);
      normalize(n_vector);
    }
  }

  void Polygon::make_uv_vectors()
  {
    for (unsigned int i = 0; i < p_3D.size(); ++i)
    {
      Vector<Real_t> dp = p_3D[i] - p_3D[0];
      Real_t u_i = u_vector * dp;
      Real_t v_i = v_vector * dp;
      u_list.push_back(u_i);
      v_list.push_back(v_i);
    }
  }

  void Polygon::make_transform_matrix()
  {
    Real_t m[3][3];
    m[0][0] = u_vector.x;
    m[0][1] = u_vector.y;
    m[0][2] = u_vector.z;
    m[1][0] = v_vector.x;
    m[1][1] = v_vector.y;
    m[1][2] = v_vector.z;
    m[2][0] = n_vector.x;
    m[2][1] = n_vector.y;
    m[2][2] = n_vector.z;

    Real_t det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                 m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                 m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    Real_t invdet = 1. / det;

    mat_inv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
    mat_inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
    mat_inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
    mat_inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
    mat_inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
    mat_inv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
    mat_inv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
    mat_inv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
    mat_inv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
  }

  void Polygon::uv_to_xyz(Real_t u_i, Real_t v_i, Vector<Real_t> &p_i)
  {
    p_i.x = mat_inv[0][0] * u_i + mat_inv[0][1] * v_i + p_3D[0].x;
    p_i.y = mat_inv[1][0] * u_i + mat_inv[1][1] * v_i + p_3D[0].y;
    p_i.z = mat_inv[2][0] * u_i + mat_inv[2][1] * v_i + p_3D[0].z;
  }

  void Polygon::xyz_to_uv(const Vector<Real_t> &p_i, Real_t &u_i, Real_t &v_i)
  {
    Vector<Real_t> dp = p_i - p_3D[0];
    u_i = u_vector * dp;
    v_i = v_vector * dp;
  }

  bool Polygon::is_inside(double u_i, double v_i) // using Jordan curve theorem,
                                                  // checks whether the point is inside the curve or not;
                                                  // it has a bug when (v_list[i] == v_list[j])
  {
    bool c = false;
    int nvert = p_3D.size();
    //     for (int i = 0; int j = nvert-1; i < nvert; j = i++)
    int j = nvert - 1;
    for (int i = 0; i < nvert; ++i)
    {
      if (((v_list[i] > v_i) != (v_list[j] > v_i)) &&
          (u_i < (u_list[j] - u_list[i]) * (v_i - v_list[i]) / (v_list[j] - v_list[i]) + u_list[i]))
        c = !c;
      j = i;
    }
    return c;
  }

  bool Polygon::in_contact(const Vector<double> &v, const double r, Vector<double> &contact_vector)
  {
    std::string s = "incomplete function:";
    s += __FILE__ + std::to_string(__LINE__) + __func__;
    output->warning(s);
    std::cout << "  " << v << r << contact_vector << std::endl;
    return false;
  }

} // shape

CAVIAR_NAMESPACE_CLOSE
