
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Implementation of M-Shake algorithm is done by Ashkan Shahmoradi
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

#ifndef CAVIAR_OBJECTS_CONSTRAINT_MSHAKE_H
#define CAVIAR_OBJECTS_CONSTRAINT_MSHAKE_H

#include "caviar/objects/constraint.h"

#include <vector>

namespace caviar {
namespace objects {

class Domain;
namespace constraint {



/**
 * This class fixes atomic bonds using M-Shake method
 * 
 * 
 */
class M_shake : public Constraint {
 public:
  M_shake (class CAVIAR *);
   ~M_shake ( );
  bool read (class caviar::interpreter::Parser *);

  void step_part_I (int);
  void step_part_II (int);
  void step_part_III (int);

  void bond_fix ();

  void verify_settings();

  inline double reset_distance (double );

  static inline int delta(int a,int b) {
	  if(a==b)return 1;
  	else return 0;
  }

  class objects::Domain *domain;

  double dt;
  double error_tolerance;
  bool initialized;

};

} //constraint
} //objects
} // namespace caviar

#endif
