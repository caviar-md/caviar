
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
      else if (string_cmp(t, "all_atom_test"))
      {
        all_atom_test = true;
      }
      else if (string_cmp(t, "rebuild_test"))
      {
        rebuild_test = true;
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

    local_cutoff = cutoff;
  }

  void Verlet_list::build_neighlist()
  {
    if (all_atom_test)
    {
      all_atom_test_function(0);
      return;
    }

    calculate_cutoff_extra();

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
    ghost_pos_old = pos_ghost;
#ifdef CAVIAR_WITH_MPI
    mpi_rank_old = atom_data->atom_struct_owned.mpi_rank;
#endif
  }

} // neighborlist

CAVIAR_NAMESPACE_CLOSE
