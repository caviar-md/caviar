
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

CAVIAR_NAMESPACE_OPEN

namespace neighborlist
{

  /**
   * A verlet list.
   */
  class Verlet_list : public Neighborlist
  {
  public:
    Verlet_list(class CAVIAR *);
    bool read(class caviar::interpreter::Parser *);
    void init();
    bool rebuild_neighlist();
    void build_neighlist();
    double dt, cutoff_extra;
    double cutoff_extra_coef;

  private:
  /**
   * position of the particles at the previous verlet list generation step
  */
  std::vector<Vector<double>> pos_old;
  std::vector<int> mpi_rank_old;

  };

} // neighborlist

CAVIAR_NAMESPACE_CLOSE

#endif
