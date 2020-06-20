
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

#ifndef CAVIAR_OBJECTS_NEIGHBORLIST_VERLETLIST_H
#define CAVIAR_OBJECTS_NEIGHBORLIST_VERLETLIST_H

#include "caviar/objects/neighborlist.h"
//#include "caviar/utility/python_utils_dec.h"

namespace caviar {
namespace objects {
namespace neighborlist {

/**
 * A verlet list.
 */
class Verlet_list : public Neighborlist {
 public:
  Verlet_list (class CAVIAR *);
  bool read (class caviar::interpreter::Parser *);
  void init ();
  bool rebuild_neighlist ();
  void build_neighlist ();
  double dt, cutoff_extra;
  double cutoff_extra_coef;
//  double cutoff; // Defined in the base class
 public:
};

//void export_py_Verlet_list ();

} //neighborlist
} //objects
} // namespace caviar

#endif
