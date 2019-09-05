
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

#ifndef CAVIAR_UTILITY_H
#define CAVIAR_UTILITY_H

#include "caviar_config.h"
#include <algorithm>
#include <vector>

/**
 * This file contains some small and useful functions used in the CAVIAR.
 */

namespace caviar {

/**
 * A function that compares two strings. It can be defined case sensitive or not.
 * Object names won't be compared by this function, so they are case sensitive.
 */
template <typename T1, typename T2>
bool string_cmp(const T1 a, const T2 b) {
#ifdef CAVIAR_SCRIPT_COMMAND_CASE_INSENSITIVE
  T1 a_lowercase = a;
  T1 b_lowercase = b;
  std::transform(a_lowercase.begin(), a_lowercase.end(), a_lowercase.begin(), ::tolower);  
  std::transform(b_lowercase.begin(), b_lowercase.end(), b_lowercase.begin(), ::tolower);  
  return (a_lowercase == b_lowercase);
#else
  return (a==b);
#endif
}

/**
 * default case insensitive string compare template function
 */
template <typename T1, typename T2>
bool string_cmp_i(const T1 a, const T2 b) {
  T1 a_lowercase = a;
  T1 b_lowercase = b;
  std::transform(a_lowercase.begin(), a_lowercase.end(), a_lowercase.begin(), ::tolower);  
  std::transform(b_lowercase.begin(), b_lowercase.end(), b_lowercase.begin(), ::tolower);  
  return (a_lowercase == b_lowercase);
}

/**
 * minimum function (returns in compile time if possible)
 */
template <typename T>
constexpr T min (T a, T b) {
  return a < b ? a : b;
}

/**
 * maximum function (returns in compile time if possible)
 */
template <typename T>
constexpr T max (T a, T b) {
  return a > b ? a : b;
}

/**
 * integer power function (returns in compile time if possible)
 */
template <typename T>
constexpr T ipow (T num, unsigned pow) {
  return pow ? num*ipow(num, pow-1) : 1;
}

/**
 *  a fast matrix inverse template function. The input-output matrices should be of the
 *  type 'std::vector < std::vector < typename > >' . It
 */
template <typename T>
int matrix_inverse(std::vector<std::vector<T>>&A, std::vector<std::vector<T>>&A_inv) {

  const auto A_size = A.size();
  for (auto &&i : A) {
    if (i.size() != A_size)
      return 3;
  }

  switch (A_size) {

  case(2) : {
    T det=(A[0][0]*A[1][1]) - (A[0][1]*A[1][0]);

    if (det==static_cast<T>(0)) {
      return 1;
    }

    T det_inv = 1.0 / det;

    A_inv[0][0] =  A[1][1] * det_inv;
    A_inv[0][1] = -A[0][1] * det_inv;
    A_inv[1][0] = -A[1][0] * det_inv;
    A_inv[1][1] =  A[0][0] * det_inv;
    return 0;
  }

  case(3) : {
    T det=(A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2])) - 
  	    (A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])) + 
	      (A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]));//+0.00001;

    if (det==static_cast<T>(0)) {
      return 1;
    }


    double det_inv = 1.0 / det;

    A_inv[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) * det_inv;
    A_inv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * det_inv;
    A_inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * det_inv;
    A_inv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * det_inv;
    A_inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * det_inv;
    A_inv[1][2] = (A[1][0] * A[0][2] - A[0][0] * A[1][2]) * det_inv;
    A_inv[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) * det_inv;
    A_inv[2][1] = (A[2][0] * A[0][1] - A[0][0] * A[2][1]) * det_inv;
    A_inv[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) * det_inv;
    return 0;
  }

//  case(4) : {  }


  default: {
    return 2;
  }

  }
}


} // namespace caviar

#endif
