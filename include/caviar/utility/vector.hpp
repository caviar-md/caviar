
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

#ifndef CAVIAR_VECTOR_H
#define CAVIAR_VECTOR_H

#include "caviar_config.h"

#include <iostream>
#include <cmath>

CAVIAR_NAMESPACE_OPEN

template <typename T>
struct Vector
{
  Vector()
  {
    x = 0.0;
    y = 0.0;
    z = 0.0;
  }
  Vector(T x_, T y_, T z_)
  {
    x = x_;
    y = y_;
    z = z_;
  }
  T x, y, z;
};

template <typename T>
Vector<T> operator+(const Vector<T> &lhs)
{
  return lhs;
}

template <typename T>
Vector<T> operator-(const Vector<T> &lhs)
{
  return Vector<T>{-lhs.x, -lhs.y, -lhs.z};
}

template <typename T>
Vector<T> operator+(const Vector<T> &lhs, const Vector<T> &rhs)
{
  return Vector<T>{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

template <typename T>
Vector<T> operator-(const Vector<T> &lhs, const Vector<T> &rhs)
{
  return Vector<T>{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

template <typename T>
T operator*(const Vector<T> &lhs, const Vector<T> &rhs)
{
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

template <typename T>
Vector<T> operator*(const Vector<T> &lhs, const double &rhs)
{
  return Vector<T>{lhs.x * rhs, lhs.y * rhs, lhs.z * rhs};
}

template <typename T>
Vector<T> operator/(const Vector<T> &lhs, const double &rhs)
{
  return Vector<T>{lhs.x / rhs, lhs.y / rhs, lhs.z / rhs};
}

template <typename T1, typename T2>
Vector<T1> operator/(const Vector<T1> &lhs, const Vector<T2> &rhs)
{
  return Vector<T1>{lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z};
}

template <typename T>
Vector<T> operator*(const double &lhs, const Vector<T> &rhs)
{
  return Vector<T>{rhs.x * lhs, rhs.y * lhs, rhs.z * lhs};
}

template <typename T>
bool operator==(const Vector<T> &lhs, const Vector<T> &rhs)
{
  return rhs.x == lhs.x && rhs.y == lhs.y && rhs.z == lhs.z;
}

template <typename T>
Vector<T> &operator+=(Vector<T> &lhs, const Vector<T> &rhs)
{
  lhs.x += rhs.x;
  lhs.y += rhs.y;
  lhs.z += rhs.z;
  return lhs;
}

template <typename T>
Vector<T> &operator-=(Vector<T> &lhs, const Vector<T> &rhs)
{
  lhs.x -= rhs.x;
  lhs.y -= rhs.y;
  lhs.z -= rhs.z;
  return lhs;
}

template <typename T>
Vector<T> &operator*=(Vector<T> &lhs, const double &rhs)
{
  lhs.x *= rhs;
  lhs.y *= rhs;
  lhs.z *= rhs;
  return lhs;
}

template <typename T>
Vector<T> &operator/=(Vector<T> &lhs, const double &rhs)
{
  lhs.x /= rhs;
  lhs.y /= rhs;
  lhs.z /= rhs;
  return lhs;
}

template <typename T>
std::ostream &operator<<(std::ostream &out, const Vector<T> &rhs)
{
  return out << rhs.x << ' ' << rhs.y << ' ' << rhs.z;
}

template <typename T>
std::istream &operator<<(std::istream &in, Vector<T> &rhs)
{
  return in >> rhs.x >> rhs.y >> rhs.z;
}

template <typename T>
constexpr Vector<T> cross_product(const Vector<T> &v1, const Vector<T> &v2)
{
  return Vector<T>{v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x};
}

template <typename T>
constexpr T dot_product(const Vector<T> &v1, const Vector<T> &v2)
{
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
constexpr T norm(const Vector<T> &v1)
{
  return std::sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
}

// https://software.intel.com/en-us/node/583439
// page 180 of: http://www.openmp.org/mp-documents/OpenMP4.0.0.pdf
#ifdef CAVIAR_WITH_OPENMP

#pragma omp declare reduction(+                          \
                              : caviar::Vector <int64_t> \
                              : omp_out += omp_in)

#pragma omp declare reduction(+                           \
                              : caviar::Vector <uint64_t> \
                              : omp_out += omp_in)

#pragma omp declare reduction(+                          \
                              : caviar::Vector <int32_t> \
                              : omp_out += omp_in)

#pragma omp declare reduction(+                         \
                              : caviar::Vector <double> \
                              : omp_out += omp_in)

#pragma omp declare reduction(+                        \
                              : caviar::Vector <float> \
                              : omp_out += omp_in)

#pragma omp declare reduction(+                           \
                              : caviar::Vector <uint32_t> \
                              : omp_out += omp_in)

#endif

CAVIAR_NAMESPACE_CLOSE

#endif
