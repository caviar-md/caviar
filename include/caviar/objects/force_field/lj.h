
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_LJ_H
#define CAVIAR_OBJECTS_FORCEFIELD_LJ_H

#include "caviar/objects/force_field.h"

namespace caviar {

namespace force_field {

/**
 * This class calculates LJ potential for the particles.
 * 
 */
class Lj : public Force_field {
public:
  Lj (class CAVIAR *);
  ~Lj () {};
  
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:
  bool input_by_array;
  std::vector<std::vector<Real_t>> epsilon,sigma;
  bool make_off_diagonal_vectors;

  bool input_by_atom;
  // epsilon - sigma of a single type. inter-type values will be deduced using these
  std::vector<Real_t> epsilon_atom, sigma_atom; 

  // some helper variables, used for debugging
  bool jump_fix, monitor_jump;
  double jump_tol;

  bool wca; //Week-Chandler-Anderson (WCA) potential activated.
  bool cutoff_list_activated;
  std::vector<std::vector<Real_t>> cutoff_list; // list of cutoffs when it is needed.
                                                // for example in WCA potentials
  
};

} //force_field

} // namespace caviar

#endif
