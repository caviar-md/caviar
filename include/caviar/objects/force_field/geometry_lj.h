
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_GEOMETRYLJ_H
#define CAVIAR_OBJECTS_FORCEFIELD_GEOMETRYLJ_H

#include "caviar/objects/force_field.h"


CAVIAR_NAMESPACE_OPEN

class Shape; 
namespace unique
{
  class Time_function_3d;
}
namespace force_field {

/**
 * This class makes a LJ force-field for geometry shapes
 *
 */
class Geometry_lj : public Force_field {
public:
  Geometry_lj (class CAVIAR *);
  ~Geometry_lj ();
  
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:  
  std::vector<caviar::Shape *> shape;

  unique::Time_function_3d *position_offset = nullptr;

  // contains the type of epsilon_wall-sigma_wall. This is the case when we have
  // more than one shape objects but we want to have the same parameters.
  std::vector<int> shape_type;

  // epsilon - sigma of a single type. inter-type values will be deduced using these
  std::vector<Real_t> epsilon_atom, sigma_atom; 
  std::vector<Real_t> epsilon_wall, sigma_wall; 

  // the epsilon-sigma of a LJ potential. first number is the shape index.
  // sigma[2][3] is the sigma of shape_2 and atom_type_3
  std::vector<std::vector<Real_t>> epsilon, sigma; 

  // a force_coef in case we need a total shift of the forces.
  Real_t force_coef;

  bool wca; //Week-Chandler-Anderson (WCA) potential activated.
  bool cutoff_list_activated;
  std::vector<std::vector<Real_t>> cutoff_list; // list of cutoffs when it is needed.
                                                // for example in WCA potentials
};

} //force_field

CAVIAR_NAMESPACE_CLOSE

#endif
