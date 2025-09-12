
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

#ifndef CAVIAR_OBJECTS_SHAPE_H
#define CAVIAR_OBJECTS_SHAPE_H

#include "caviar/utility/objects_common_headers.h"

CAVIAR_NAMESPACE_OPEN

inline void normalize(Vector<Real_t> &v)
{
  v /= std::sqrt(v * v);
}

/**
 * This class is the base class for all the shapes.
 *
 *
 */
class Shape : public Pointers
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
  virtual bool is_inside(const Vector<double> &) = 0;
  virtual bool is_outside(const Vector<double> &);
  virtual bool is_inside(const Vector<double> &, const double rad) = 0;
  virtual bool is_outside(const Vector<double> &, const double rad);
  virtual bool in_contact(const Vector<double> &, const double rad, Vector<double> &contact_vector) = 0;

  /**
   * Used in barostat scaling for geometrical forces.
   * 
  */   
  virtual void scale_position(double scale_ratio, caviar::Vector<int> scale_axis);
  
  FC_BASE_OBJECT_COMMON_TOOLS
};

CAVIAR_NAMESPACE_CLOSE

#endif
