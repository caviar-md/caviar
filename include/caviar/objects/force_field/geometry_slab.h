
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_GEOMETRYSLAB_H
#define CAVIAR_OBJECTS_FORCEFIELD_GEOMETRYSLAB_H

#include "caviar/objects/force_field.h"

namespace caviar {

namespace unique
{
  class Time_function_3d;
}
namespace force_field {

/**
 * This class creates a force-field for the slab geometries.
 * The default is symmetric, meaning that the slab has no direction
 * The assymetric slab means that the position of the particles relative to 
 * the slab direction matters. It is good for using soft forcefields.
 */
class Geometry_slab : public Force_field {
public:
  Geometry_slab (class CAVIAR *);
  ~Geometry_slab ();
  
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:  
  unique::Time_function_3d *position_offset = nullptr;
  unique::Time_function_3d *velocity_offset = nullptr;
  int symmetric;
  int slab_direction;
  double slab_position;
  std::vector<double> radius;
  double young_modulus, dissip_coef;
};

} //force_field

} // namespace caviar

#endif
