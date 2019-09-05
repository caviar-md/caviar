
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

#ifndef CAVIAR_OBJECTS_NEIGHBORLIST_H
#define CAVIAR_OBJECTS_NEIGHBORLIST_H

#include "caviar/utility/objects_common_headers.h"

namespace caviar {

namespace objects {
class Atom_data;

/**
 * This class is the base class for all the neighborlists.
 * Neighborlists are a way to reduce the order of MD force calculations
 * 
 */ 
class Neighborlist : public Pointers {

public:
  Neighborlist (class CAVIAR *);
  virtual ~Neighborlist ();
  virtual bool read (class caviar::interpreter::Parser *) = 0;
  virtual void init () = 0;
  virtual bool rebuild_neighlist () = 0;
  virtual void build_neighlist () = 0;
  std::vector<std::vector<LocalID_t>> neighlist;

// 'Cell_list' public functions and variables;
  virtual Vector<int> binlist_index (const Vector<double> &);
  virtual int neigh_bin_index (const Vector<double> &);
  std::vector<std::vector<std::vector<std::vector<int>>>> binlist;
  std::vector<std::vector<Vector<int>>> neigh_bin;
  Vector<int> no_bins;
  double cutoff; 

  class caviar::objects::Atom_data *atom_data;

  FC_BASE_OBJECT_COMMON_TOOLS
};

} //objects

} // namespace caviar

#endif
