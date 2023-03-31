
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_GEOMETRYSPHERELJ_H
#define CAVIAR_OBJECTS_FORCEFIELD_GEOMETRYSPHERELJ_H

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
class Geometry_sphere_lj : public Force_field {
public:
  Geometry_sphere_lj (class CAVIAR *);
  ~Geometry_sphere_lj ();
  
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:  
  unique::Time_function_3d *position_offset = nullptr;

  bool inside;
  double radius;
  caviar::Vector<double> center;

  // epsilon - sigma of a single type. inter-type values will be deduced using these
  std::vector<Real_t> epsilon_atom, sigma_atom; 
  double epsilon_wall, sigma_wall; 

  // the epsilon-sigma of a LJ potential.
  std::vector<Real_t> epsilon, sigma; 

  // a force_coef in case we need a total shift of the forces.
  Real_t force_coef;

  bool wca; //Week-Chandler-Anderson (WCA) potential activated.
  bool cutoff_list_activated;
  std::vector<Real_t> cutoff_list; // list of cutoffs when it is needed.
                                                // for example in WCA potentials
  
  
};

} //force_field

} // namespace caviar

#endif
