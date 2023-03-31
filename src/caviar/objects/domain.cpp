
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

#include "caviar/objects/domain.h"
#include "caviar/interpreter/communicator.h"

CAVIAR_NAMESPACE_OPEN

Domain::Domain(CAVIAR *fptr) : Pointers{fptr},
                               boundary_condition{Vector<int>{0, 0, 0}},
                               grid_index_x{0}, grid_index_y{0}, grid_index_z{0},
                               nprocs_x{1}, nprocs_y{1}, nprocs_z{1}, me{0}, nprocs{1}
{
  FC_OBJECT_INITIALIZE
#if defined(CAVIAR_WITH_MPI)
  MPI_Comm_rank(mpi_comm, &me);
  MPI_Comm_size(mpi_comm, &nprocs);
#endif
}

Domain::~Domain()
{
}

void Domain::verify_settings()
{
}

void Domain::calculate_local_domain()
{
#if defined(CAVIAR_WITH_MPI)
  lower_local.x = lower_global.x + (upper_global.x - lower_global.x) * grid_index_x / nprocs_x;
  lower_local.y = lower_global.y + (upper_global.y - lower_global.y) * grid_index_y / nprocs_y;
  lower_local.z = lower_global.z + (upper_global.z - lower_global.z) * grid_index_z / nprocs_z;

  upper_local.x = lower_global.x + (upper_global.x - lower_global.x) * (grid_index_x + 1) / nprocs_x;
  upper_local.y = lower_global.y + (upper_global.y - lower_global.y) * (grid_index_y + 1) / nprocs_y;
  upper_local.z = lower_global.z + (upper_global.z - lower_global.z) * (grid_index_z + 1) / nprocs_z;
#endif
}

int Domain::grid2rank(int x, int y, int z)
{
  return x + y * nprocs_x + z * nprocs_x * nprocs_y;
} // calculates process rank from it grid index

void Domain::calculate_procs_grid()
{
#ifdef CAVIAR_WITH_MPI
  neighborlist_domains.push_back(me); // don't know if it's neccesary
  find_best_grid();
  if (me == 0)
    std::cout << "info: Number of processes:\n\tnpx:" << nprocs_x << " npy:"
              << nprocs_y << " npz:" << nprocs_z << std::endl;

  grid_index_x = me % nprocs_x;
  grid_index_y = me / nprocs_x % nprocs_y;
  grid_index_z = me / nprocs_x / nprocs_y;

  for (auto i = -1; i < 2; ++i)
    for (auto j = -1; j < 2; ++j)
      for (auto k = -1; k < 2; ++k)
      {
        int ii = 0, jj = 0, kk = 0; // Zero value is just to stop '-Wmaybe-uninitialized'
        if (i == 0)
          ii = grid_index_x;
        if (i == -1)
          ii = grid_index_x == 0 ? nprocs_x - 1 : (grid_index_x - 1) % nprocs_x;
        if (i == +1)
          ii = grid_index_x == nprocs_x - 1 ? 0 : (grid_index_x + 1) % nprocs_x;
        if (j == 0)
          jj = grid_index_y;
        if (j == -1)
          jj = grid_index_y == 0 ? nprocs_y - 1 : (grid_index_y - 1) % nprocs_y;
        if (j == +1)
          jj = grid_index_y == nprocs_y - 1 ? 0 : (grid_index_y + 1) % nprocs_y;
        if (k == 0)
          kk = grid_index_z;
        if (k == -1)
          kk = grid_index_z == 0 ? nprocs_z - 1 : (grid_index_z - 1) % nprocs_z;
        if (k == +1)
          kk = grid_index_z == nprocs_z - 1 ? 0 : (grid_index_z + 1) % nprocs_z;
        all[i + 1][j + 1][k + 1] = grid2rank(ii, jj, kk);

        bool domain_found = false;
        for (auto l : neighborlist_domains)
          if (all[i + 1][j + 1][k + 1] == l)
            domain_found = true;
        if (!domain_found)
          neighborlist_domains.push_back(all[i + 1][j + 1][k + 1]);
      }
#else
  neighborlist_domains.push_back(0); // don't know if it's neccesary
#endif
}

#if defined(CAVIAR_WITH_MPI)
static std::vector<std::array<int, 3>> possible_grids(int nprocs)
{
  std::vector<std::array<int, 3>> grids;
  for (auto nprocs_x = 1; nprocs_x <= nprocs; ++nprocs_x)
  {
    if (nprocs % nprocs_x)
      continue;
    auto nprocs_yz = nprocs / nprocs_x;
    for (auto nprocs_y = 1; nprocs_y <= nprocs_yz; ++nprocs_y)
    {
      if (nprocs_yz % nprocs_y)
        continue;
      grids.push_back({nprocs_x, nprocs_y, nprocs_yz / nprocs_y});
    }
  }
  return grids;
}
#endif

void Domain::find_best_grid()
{
#if defined(CAVIAR_WITH_MPI)
  auto box_length_x = upper_global.x - lower_global.x;
  auto box_length_y = upper_global.y - lower_global.y;
  auto box_length_z = upper_global.z - lower_global.z;

  auto area_xy = box_length_x * box_length_y, area_xz = box_length_x * box_length_z, area_yz = box_length_y * box_length_z;

  auto grids = possible_grids(nprocs);
  auto min_area = 1.1 * (area_xy + area_xz + area_yz);
  auto min_area_index = 0;
  for (unsigned int i = 0; i < grids.size(); ++i)
  {
    auto area = area_xy / grids[i][0] / grids[i][1] + area_xz / grids[i][0] / grids[i][2] + area_yz / grids[i][1] / grids[i][2];
    if (area < min_area)
    {
      min_area = area;
      min_area_index = i;
    }
  }
  nprocs_x = grids[min_area_index][0];
  nprocs_y = grids[min_area_index][1];
  nprocs_z = grids[min_area_index][2];
#endif
}

// note that this force cannot be used when one is using neighborlist loop.
// currently its main usage is when one is using Shake like algorithms.
// or Spring_bond or Spring_angle force_fields.
Vector<Real_t> Domain::periodic_distance(const Vector<Real_t> v)
{
  caviar::Vector<Real_t> vf = v;
  static caviar::Vector<double> domain_dh = {0.5 * (upper_global.x - lower_global.x),
                                             0.5 * (upper_global.y - lower_global.y),
                                             0.5 * (upper_global.z - lower_global.z)};
  if (boundary_condition.x == 1)
  {
    while (vf.x > +domain_dh.x)
    {
      vf.x -= domain_dh.x * 2.0;
    }
    while (vf.x < -domain_dh.x)
    {
      vf.x += domain_dh.x * 2.0;
    }
  }

  if (boundary_condition.y == 1)
  {
    while (vf.y > +domain_dh.y)
    {
      vf.y -= domain_dh.y * 2.0;
    }
    while (vf.y < -domain_dh.y)
    {
      vf.y += domain_dh.y * 2.0;
    }
  }

  if (boundary_condition.z == 1)
  {
    while (vf.z > +domain_dh.z)
    {
      vf.z -= domain_dh.z * 2.0;
    }
    while (vf.z < -domain_dh.z)
    {
      vf.z += domain_dh.z * 2.0;
    }
  }

  // std::cout << domain_dh << std::endl;
  // if( !( v.x == vf.x && v.y == vf.y && v.z == vf.z) )
  // std::cout << "x: " << v << " vs " << vf << std::endl;
  return vf;
}

CAVIAR_NAMESPACE_CLOSE
