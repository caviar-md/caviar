
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
 * A cell list with an internal verlet list.
 * to this point, the
 * cell list usage is in the potential calculation of electrostatic_ewald_r force_field
 * on the finite_element mesh boundaries,
 */
class Neighborlist : public Pointers
{
public:
  Neighborlist(class CAVIAR *);

  virtual ~Neighborlist();

  bool read(class caviar::interpreter::Parser *);

  virtual void init();

  /**
   *  if update_neighborlist == true, it makes a new list,
   *  if not, it checks if a new list is needed, and if the 
   *  result is true, it creates a new list.
   *  It returns weather or not a new list is created
   */
  virtual bool build(bool update_neighborlist);

  /**
   * Check if the update is needed due to particle movements
  */
  virtual bool update_is_needed();

  virtual void update_verlet_list();

  virtual void update_verlet_list_from_cell_list();

  virtual void update_cell_list();

  virtual void calculate_cutoff_extra();

  virtual void all_atom_test_function(int state = 0);

  virtual Vector<int> binlist_index(const Vector<double> &);

  virtual int neigh_bin_index(const Vector<double> &);

  virtual void make_neigh_bin();

  class Domain *domain = nullptr;

  class caviar::Atom_data *atom_data = nullptr;

  /**
   * if 'true' one can use this class as a 'verlet_list'.
   */
  bool make_verlet_list_from_cell_list = false;

  /**
   *
   * Verlet list: the list of the neighbors of each atom
   */
  std::vector<std::vector<unsigned int>> neighlist;

  /**
   * Cell List: The list of particles index in each bin with 3 Index in 3D
   */
  std::vector<std::vector<std::vector<std::vector<unsigned int>>>> binlist;

  ///**
  // * Cell List: The list of particles index in each bin
  // */
  //std::vector<std::vector<unsigned int>> binlist_linear;


  /**
   * Cell List: The 3D index of the neighbors of each bin.
   */
  std::vector<std::vector<Vector<int>>> neigh_bin;

  ///*
  // * Cell List: The 3D index of the neighbors of each bin.
  // */
  //std::vector<std::vector<unsigned int>> neigh_bin_linear;

  /**
   * Cell List: Number of bins in each direction
   */
  Vector<int> no_bins;

  /**
   * Maximum cutoff of short-ranged force-fields.
   */
  double cutoff = 0;

  /**
   * cutoff_extra = cutoff_extra_coef * cutoff
   */
  double cutoff_extra = 0;

  /**
   * square of rebuild verlet list threshold distance ~ (0.5*(cutoff_extra - cutoff))^2
   */
  double threshold_distance_sq;

  /**
   * default extra coef
   */
  double cutoff_extra_coef = 1.12347;

  /**
   * Timestep, Used in cutoff_extra calculations
   */
  double dt;

  /**
   * MPI rank of the classs
   */
  int my_mpi_rank = -1;

protected:
  /**
   *  default type is verlet_list. Celllist is built if reqired by some force_field
   *  since it is more computational intensive thatn verlet list
   */
  double build_cell_list = false;

  /**
   *  inverse of cutoff. For faster computations
   */
  double cutoff_inv;

  /**
   * position of the ghost particles at the previous verlet list generation step
   */
  std::vector<Vector<double>> ghost_pos_old;

  /**
   * position of the particles at the previous verlet list generation step
   */
  std::vector<Vector<double>> pos_old;

  std::vector<int> mpi_rank_old;

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

  bool initialized = false;

public:
  FC_BASE_OBJECT_COMMON_TOOLS
};

CAVIAR_NAMESPACE_CLOSE

#endif
