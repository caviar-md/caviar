
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_GEOMETRY_H
#define CAVIAR_OBJECTS_FORCEFIELD_GEOMETRY_H

#include "caviar/objects/force_field.h"

CAVIAR_NAMESPACE_OPEN

class Shape;
namespace unique
{
  class Time_function_3d;
}
namespace force_field
{

  /**
   * This class creates a spring force-field for the shape geometries
   *
   */
  class Geometry : public Force_field
  {
  public:
    Geometry(class CAVIAR *);
    ~Geometry();

    bool read(class caviar::interpreter::Parser *);
    void verify_settings();
    void calculate_acceleration();

  public:
    unique::Time_function_3d *position_offset = nullptr;
    unique::Time_function_3d *velocity_offset = nullptr;
    std::vector<caviar::Shape *> shape;
    bool shape_size_warning;
    std::vector<double> radius;
    double young_modulus, dissip_coef;
  };

} // force_field

CAVIAR_NAMESPACE_CLOSE

#endif
