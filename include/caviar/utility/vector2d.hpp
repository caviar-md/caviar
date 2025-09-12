
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



#include <iostream>
#include <cmath>

namespace caviar
{

  template <typename T>
  struct Vector2d
  {
    Vector2d()
    {
      x = 0.0;
      y = 0.0;
    }
    Vector2d(T x_, T y_)
    {
      x = x_;
      y = y_;
    }
    Vector2d(const Vector2d<T> &v)
    {
      x = v.x;
      y = v.y;
    }
    T x, y;
  };

  template <typename T>
  Vector2d<T> operator+(const Vector2d<T> &lhs)
  {
    return lhs;
  }

  template <typename T>
  Vector2d<T> operator-(const Vector2d<T> &lhs)
  {
    return Vector2d<T>{-lhs.x, -lhs.y};
  }

  template <typename T>
  Vector2d<T> operator+(const Vector2d<T> &lhs, const Vector2d<T> &rhs)
  {
    return Vector2d<T>{lhs.x + rhs.x, lhs.y + rhs.y};
  }

  template <typename T>
  Vector2d<T> operator-(const Vector2d<T> &lhs, const Vector2d<T> &rhs)
  {
    return Vector2d<T>{lhs.x - rhs.x, lhs.y - rhs.y};
  }

  template <typename T>
  T operator*(const Vector2d<T> &lhs, const Vector2d<T> &rhs)
  {
    return lhs.x * rhs.x + lhs.y * rhs.y;
  }

  template <typename T>
  Vector2d<T> operator*(const Vector2d<T> &lhs, const double &rhs)
  {
    return Vector2d<T>{lhs.x * rhs, lhs.y * rhs};
  }

  template <typename T>
  Vector2d<T> operator/(const Vector2d<T> &lhs, const double &rhs)
  {
    return Vector2d<T>{lhs.x / rhs, lhs.y / rhs};
  }

  template <typename T>
  Vector2d<T> operator*(const double &lhs, const Vector2d<T> &rhs)
  {
    return Vector2d<T>{rhs.x * lhs, rhs.y * lhs};
  }

  template <typename T>
  bool operator==(const Vector2d<T> &lhs, const Vector2d<T> &rhs)
  {
    return rhs.x == lhs.x && rhs.y == lhs.y;
  }

  template <typename T>
  Vector2d<T> &operator+=(Vector2d<T> &lhs, const Vector2d<T> &rhs)
  {
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    return lhs;
  }

  template <typename T>
  Vector2d<T> &operator-=(Vector2d<T> &lhs, const Vector2d<T> &rhs)
  {
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    return lhs;
  }

  template <typename T>
  Vector2d<T> &operator*=(Vector2d<T> &lhs, const double &rhs)
  {
    lhs.x *= rhs;
    lhs.y *= rhs;
    return lhs;
  }

  template <typename T>
  Vector2d<T> &operator/=(Vector2d<T> &lhs, const double &rhs)
  {
    lhs.x /= rhs;
    lhs.y /= rhs;
    return lhs;
  }

  template <typename T>
  std::ostream &operator<<(std::ostream &out, const Vector2d<T> &rhs)
  {
    return out << rhs.x << '\t' << rhs.y;
  }

  template <typename T>
  std::istream &operator<<(std::istream &in, Vector2d<T> &rhs)
  {
    return in >> rhs.x >> rhs.y;
  }

  template <typename T>
  constexpr T dot_product(const Vector2d<T> &v1, const Vector2d<T> &v2)
  {
    return v1.x * v2.x + v1.y * v2.y;
  }

  template <typename T>
  constexpr T norm(const Vector2d<T> &v1)
  {
    return std::sqrt(v1.x * v1.x + v1.y * v1.y);
  }

}
