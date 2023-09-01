
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

#include "caviar/objects/neighborlist/verlet_list.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"

#include <cmath>

CAVIAR_NAMESPACE_OPEN

namespace neighborlist
{
  Verlet_list::Verlet_list(CAVIAR *fptr) : Neighborlist{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    cutoff = -1.0;
    dt = -1.0;
    cutoff_extra = 0.0; // this value changes every timestep
    cutoff_extra_coef = 0.0;
  }

  bool Verlet_list::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else if (string_cmp(t, "cutoff"))
      {
        GET_OR_CHOOSE_A_REAL(cutoff, "", "")
        if (cutoff < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "cutoff have to non-negative.");
      }
      else if (string_cmp(t, "dt"))
      {
        GET_OR_CHOOSE_A_REAL(dt, "", "")
        if (dt < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "dt have to non-negative.");
      }
      else if (string_cmp(t, "cutoff_extra_coef"))
      {
        GET_OR_CHOOSE_A_REAL(cutoff_extra_coef, "", "")
        if (cutoff_extra_coef < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "cutoff_extra_coef have to non-negative.");
      }
      else
        error->all(FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
    }

    return in_file;
  }

  void Verlet_list::init()
  {
    FC_NULLPTR_CHECK(atom_data)
    my_mpi_rank = atom_data->get_mpi_rank();

    if (cutoff <= 0.0)
      error->all(FC_FILE_LINE_FUNC, "cutoff have to non-negative.");
    if (dt <= 0.0)
      error->all(FC_FILE_LINE_FUNC, "dt have to non-negative.");

    const auto &pos = atom_data->atom_struct_owned.position;
    pos_old.resize(pos.size());
  }

  bool Verlet_list::rebuild_neighlist()
  {

    const auto &pos = atom_data->atom_struct_owned.position; 

    unsigned int pos_size = pos.size();

    if (pos_size != pos_old.size())
    {
      return true;
    }
    else
    {
      for (unsigned int i = 0; i < pos_size; ++i)
      {
#ifdef CAVIAR_WITH_MPI
        if (atom_data->atom_struct_owned.mpi_rank[i] != mpi_rank_old[i]) //make new verlet_list if any mpi_rank is changed,
          return true;

        if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank) // then, ignore mpi_ranks out of the domain
          continue;
#endif

        auto disp = pos[i] - pos_old[i];
        if (disp * disp > cutoff_extra * cutoff_extra / 4)
          return true;
      }
    }

    return false;
  }

  void Verlet_list::build_neighlist()
  {

    if (cutoff_extra_coef > 0)
    {
      const auto &vel = atom_data->atom_struct_owned.velocity;

      double max_vel_sq = 0.0;
      for (unsigned int i = 0; i < vel.size(); ++i)
      {
        double vel_sq_temp = vel[i] * vel[i];
        if (max_vel_sq < vel_sq_temp)
          max_vel_sq = vel_sq_temp;
      }
      cutoff_extra = cutoff_extra_coef * dt * std::sqrt(max_vel_sq);
    }

    const auto &pos = atom_data->atom_struct_owned.position;
    const auto &pos_ghost = atom_data->atom_struct_ghost.position;
    const auto cutoff_sq = (cutoff + cutoff_extra) * (cutoff + cutoff_extra);
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
        auto dr = pos[j] - pos[i];
        if (dr * dr < cutoff_sq)
        {
          neighlist[i].emplace_back(j);
        }
      }
      for (unsigned int j = 0; j < pos_ghost.size(); ++j)
      {
        auto dr = pos_ghost[j] - pos[i];
        if (dr * dr < cutoff_sq)
        {
          neighlist[i].emplace_back(j + pos_size);
        }
      }
    }

    pos_old = pos;
    mpi_rank_old = atom_data->atom_struct_owned.mpi_rank;
  }

} // neighborlist

CAVIAR_NAMESPACE_CLOSE
