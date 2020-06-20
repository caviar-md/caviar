
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

#ifndef CAVIAR_OBJECTS_NEIGHBORLIST_CELLLIST_H
#define CAVIAR_OBJECTS_NEIGHBORLIST_CELLLIST_H

#include "caviar/objects/neighborlist.h"
//#include "caviar/utility/python_utils_dec.h"

namespace caviar {
namespace objects {
class Domain;
namespace neighborlist {

/**
 * A cell list with an internal verlet list.
 * to this point, the
 * cell list usage is in the potential calculation of electrostatic_ewald_r force_field
 * on the finite_element mesh boundaries, 
 */
class Cell_list : public Neighborlist {
 public:
  Cell_list (class CAVIAR *);
  bool read (class caviar::interpreter::Parser *);
  void init ();
  bool rebuild_neighlist ();
  void build_neighlist ();
  void build_binlist ();
  Vector<int> binlist_index (const Vector<double> &);
  int neigh_bin_index (const Vector<double> &);
 public:
  void make_neigh_bin ();
  std::shared_ptr<class objects::Domain > domain;
  bool domain_set;

  /**
   * if 'true' one can use this class as a 'verlet_list'.
   */
  bool make_neighlist; 

  double cutoff_neighlist;
};

//void export_py_Cell_list ();

} //neighborlist
} //objects
} // namespace caviar

#endif
