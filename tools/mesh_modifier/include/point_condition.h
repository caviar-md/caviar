
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

#ifndef POINT_CONDITION_H
#define POINT_CONDITION_H

#include <string>
#include "vector.h"
class Point_condition
{
public:
  Point_condition();
  Point_condition(const std::string el, const std::string op, double);
  Point_condition(const std::string el, const Vector<double>, const std::string op, double);
  ~Point_condition();

  // to be used in functions.
  bool in_condition(const Vector<double> &) const;

private:
  // setting the operator into the given reference int variable
  void set_operator(const std::string, int &type);

  // -2: smaller, -1: equal-smaller, 0: equal, 1: equal-larger, 2: larger
  int operation_type;

  // 0: x, 1: y, 2: z
  int vector_component;

  // 0: compare one component 'vector_component' using 'operation_type' with 'value',
  // 1: compare the cartesian distance from the point to the 'vec'
  int condition_type;

  // the scalar value to be compared with
  double value, value_sq;

  // Vector to be used and-or compared
  Vector<double> vec;
};

#endif
