
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

#include "caviar/utility/objects_common_headers.hpp"

namespace caviar
{

  inline void normalize(Vector3d<double> &v)
  {
    v /= std::sqrt(v * v);
  }

  /**
   * This class is the base class for all the shapes.
   *
   *
   */
  class Shape 
  {
  public:
    /**
     * Constructor.
     */
    Shape(class CAVIAR *);

    /**
     * Destructor.
     */
    virtual ~Shape();

    virtual bool read(class caviar::interpreter::Parser *) = 0;
    virtual bool is_inside(const Vector3d<double> &) = 0;
    virtual bool is_outside(const Vector3d<double> &);
    virtual bool is_inside(const Vector3d<double> &, const double rad) = 0;
    virtual bool is_outside(const Vector3d<double> &, const double rad);
    virtual bool in_contact(const Vector3d<double> &, const double rad, Vector3d<double> &contact_vector) = 0;

    /**
     * Used in barostat scaling for geometrical forces.
     *
     */
    virtual void scale_position(double scale_ratio, caviar::Vector3d<int> scale_axis);

    FC_BASE_OBJECT_COMMON_TOOLS
  };

}
