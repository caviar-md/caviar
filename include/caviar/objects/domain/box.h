
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

#ifndef CAVIAR_OBJECTS_DOMAIN_BOX_H
#define CAVIAR_OBJECTS_DOMAIN_BOX_H

#include "caviar/objects/domain.h"
#include "caviar/utility/python_utils_dec.h"

namespace caviar {
namespace objects {
namespace domain {

/**
 * This class emulates a 3D box for a simulation domain
 * 
 * 
 */
class Box : public Domain {
public:
  Box (class CAVIAR *);
  bool read (class caviar::interpreter::Parser *);
  void calculate_local_domain ();
  void generate ();  

  double fix_distance_x(double d);
  double fix_distance_y(double d);
  double fix_distance_z(double d);

  Vector<double> fix_distance(caviar::Vector<double> v); 

  Vector <Real_t> half_edge;


  FC_PYDEC_SETGET_CAVVEC(half_edge,Real_t);  

  FC_PYDEC_SETGET_CAVVEC(upper_local,Real_t);  
  FC_PYDEC_SETGET_CAVVEC(lower_local,Real_t);  
  FC_PYDEC_SETGET_CAVVEC(upper_global,Real_t);  
  FC_PYDEC_SETGET_CAVVEC(lower_global,Real_t);  

  FC_PYDEC_SETGET_CAVVEC(boundary_condition,int);  

public:
  
};

void export_py_Box ();

} //domain
} //objects
} // namespace caviar

#endif
