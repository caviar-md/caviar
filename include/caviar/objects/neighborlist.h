
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

CAVIAR_NAMESPACE_OPEN

class Atom_data;

/**
 * This class is the base class for all the neighborlists.
 * Neighborlists are a way to reduce the order of MD force calculations
 *
 */
class Neighborlist : public Pointers
{

public:
  Neighborlist(class CAVIAR *);
  virtual ~Neighborlist();
  virtual bool read(class caviar::interpreter::Parser *) = 0;
  virtual void init() = 0;
  virtual bool rebuild_neighlist();
  virtual void build_neighlist() = 0;
  virtual void calculate_cutoff_extra();
  virtual void all_atom_test_function(int state=0);
  std::vector<std::vector<LocalID_t>> neighlist;

  // 'Cell_list' public functions and variables;
  virtual Vector<int> binlist_index(const Vector<double> &);
  virtual int neigh_bin_index(const Vector<double> &);
  std::vector<std::vector<std::vector<std::vector<int>>>> binlist;
  std::vector<std::vector<Vector<int>>> neigh_bin;
  Vector<int> no_bins;
  double dt;
  double cutoff;
  /**
   * if (cutoff_extra_coef > 0 ) cutoff_extra calculated according to maximum velocity of the particles and added to cutoff in order 
   * to have less neighborlist re-make.
  */
  double cutoff_extra; 
  double cutoff_extra_coef;

  class caviar::Atom_data *atom_data;
  /**
   * MPI rank of the classs
  */
  int my_mpi_rank = -1;

  protected:

  /**
   * position of the ghost particles at the previous verlet list generation step
  */
  std::vector<Vector<double>> ghost_pos_old;

  /**
   * position of the particles at the previous verlet list generation step
  */
  std::vector<Vector<double>> pos_old;
  std::vector<int> mpi_rank_old;
  double local_cutoff; // in verlet_list, it is the same as cutoff. in cell_list, it is cutoff_neighlist;
  /**
  *  if true, all of the atoms will see all of the other atoms
  *  the result of simulation must be exactly the same as the neighborlist with correct parameters.
  */
  bool all_atom_test = false;

  /**
  *  if true, verlet_list will be made on every timestep
  *  the result of simulation must be exactly the same as the neighborlist with correct parameters.
  */
  bool rebuild_test = false;

public:

  FC_BASE_OBJECT_COMMON_TOOLS
};

CAVIAR_NAMESPACE_CLOSE

#endif
