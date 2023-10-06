
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

#ifndef CAVIAR_OBJECTS_DOMAIN_H
#define CAVIAR_OBJECTS_DOMAIN_H

#include "caviar/utility/objects_common_headers.h"

CAVIAR_NAMESPACE_OPEN

/**
 * This class is the base class for all the domains.
 * Domain contains the value of the simulation boxes.
 *
 */
class Domain : public Pointers
{
public:
  Domain(class CAVIAR *);
  virtual ~Domain();

  virtual bool read(class caviar::interpreter::Parser *);

  virtual void calculate_local_domain();
  virtual void generate();
  virtual void calculate_procs_grid();

  /**
   * It must be called after the changes in the domain
  */
  virtual void update_after_domain_change();

  /**
   * Used in barostat scaling for geometrical forces.
   *
   */
  virtual void scale_position(double scale_ratio, caviar::Vector<int> scale_axis);

  /**
   * calculates process rank from it grid index
   */
  virtual int grid2rank(int x, int y, int z);

  /**
   * with respect to the number of processes, makes a grid with lowest shared area possible.
   */
  virtual void find_best_grid();

  virtual double fix_distance_x(double d);
  virtual double fix_distance_y(double d);
  virtual double fix_distance_z(double d);

  virtual caviar::Vector<double> fix_distance(caviar::Vector<double> v);

  /**
   * Gives the corrected position in periodic boundary condition. Use case: if an atom of a molecule crossed the  global boundary
   * in periodic condition, its actual location is on the other side of the boundary, but since exchanging just one atom of a molecule
   * between domain is MPI expensive, we use fix_position to get the correct location of the atoms in cases it is needed, such as
   * when using the dealii force_field and we need all the particles remain inside mesh.
   */
  virtual caviar::Vector<double> fix_position(caviar::Vector<double> v, caviar::Vector<int> &msd_value, bool &update_verlet_list);

  virtual Vector<Real_t> periodic_distance(const Vector<Real_t>);

  /**
   * Total volume of (mpi) local domain. In non-mpi simulations, local == global
   */
  virtual double volume_local();

  /**
   * Total volume of global domain
   */
  virtual double volume_global();

  Vector<int> boundary_condition;

  int grid_index_x, grid_index_y, grid_index_z; // starts from (0) to (nprocs_i-1) ; i=x,y,z
  int nprocs_x, nprocs_y, nprocs_z;             // it can be at least (1) and at most (nprocs)

  Vector<Real_t> lower_global, upper_global;
  Vector<Real_t> lower_local, upper_local;

  Vector<Real_t> size_local, size_global;

  /**
   * used in MD_MPI case:
   * MPI process rank and number of processes
   */
  int me, nprocs;

  /**
   *  all neighborlist domains around the me=all [1][1][1]. left=all[0][1][1]. up=all[1][2][1].
   *  if one domain exists, in one process case, all[i][j][k]=me for 0<=i,j,k<=2
   *  left&right: x direction, down&up: y direction, bottom&top: z direction. right&up&top are the positive directions
   */
  int all[3][3][3];

  /**
   * defined to have a faster MPI when we have less than 27 processes or when domains can have similar neigbors
   */
  std::vector<int> neighborlist_domains;

  Vector<Real_t> half_edge;

public:
  FC_BASE_OBJECT_COMMON_TOOLS
};

CAVIAR_NAMESPACE_CLOSE

#endif
