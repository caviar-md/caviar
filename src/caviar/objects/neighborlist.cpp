
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

#include "caviar/objects/neighborlist.h"
#include "caviar/interpreter/error.h"
#include "caviar/objects/atom_data.h"

CAVIAR_NAMESPACE_OPEN

Neighborlist::Neighborlist(CAVIAR *fptr) : Pointers{fptr},
                                           atom_data{nullptr} {
                                               FC_OBJECT_INITIALIZE}

                                           Neighborlist::~Neighborlist()
{
}

void Neighborlist::verify_settings()
{
}

Vector<int> Neighborlist::binlist_index(const Vector<double> &a)
{
  error->all("NOT implemented yet");
  return Vector<int>{static_cast<int>(a.x), 0, 0};
}
int Neighborlist::neigh_bin_index(const Vector<double> &a)
{
  error->all("NOT implemented yet");
  return a.x;
}

bool Neighborlist::rebuild_neighlist()
{
  if (initialize)
  {
    initialize = false;
    init();
  }

  if (rebuild_test)
  {
    return true;
  }

  const auto &pos = atom_data->atom_struct_owned.position;
  const auto &ghost_pos = atom_data->atom_struct_ghost.position;

  unsigned int pos_size = pos.size();
  unsigned int ghost_size = ghost_pos.size();

  if (pos_size != pos_old.size() || ghost_size != ghost_pos_old.size())
  {
    return true;
  }
  else
  {
    for (unsigned int i = 0; i < pos_size; ++i)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != mpi_rank_old[i]) // make new verlet_list if any mpi_rank is changed,
        return true;

      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank) // then, ignore mpi_ranks out of the domain
        continue;
#endif

      auto disp = pos[i] - pos_old[i];
      if (disp * disp > threshold_distance_sq)
        return true;
    }

    for (unsigned int i = 0; i < ghost_size; ++i)
    {
      auto disp = ghost_pos[i] - ghost_pos_old[i];
      if (disp * disp > threshold_distance_sq)
        return true;
    }
  }

  return false;
}

void Neighborlist::calculate_cutoff_extra()
{
  // if (cutoff_extra_coef <= 0)
  // {
  //   cutoff_extra = 0;
  //   return;
  // }

  // const auto &vel = atom_data->atom_struct_owned.velocity;

  // double max_vel_sq = 0.0;
  // for (unsigned int i = 0; i < vel.size(); ++i)
  // {
  //   double vel_sq_temp = vel[i] * vel[i];
  //   if (max_vel_sq < vel_sq_temp)
  //     max_vel_sq = vel_sq_temp;
  // }
  // cutoff_extra = cutoff_extra_coef * dt * std::sqrt(max_vel_sq);
}

void Neighborlist::all_atom_test_function(int state)
{
  if (state == 0)
  {
    const auto &pos = atom_data->atom_struct_owned.position;
    const auto &pos_ghost = atom_data->atom_struct_ghost.position;
    const auto pos_size = pos.size();
    neighlist.clear();
    neighlist.resize(pos_size);

    for (unsigned int i = 0; i < pos_size; ++i)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      for (unsigned int j = i + 1; j < pos_size; ++j)
      {
#ifdef CAVIAR_WITH_MPI
        if (atom_data->atom_struct_owned.mpi_rank[j] != my_mpi_rank)
          continue;
#endif
        neighlist[i].emplace_back(j);
      }
      for (unsigned int j = 0; j < pos_ghost.size(); ++j)
      {
        neighlist[i].emplace_back(j + pos_size);
      }
    }
  }
}
CAVIAR_NAMESPACE_CLOSE
