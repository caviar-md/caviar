
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

#include "caviar/objects/force_field/lj_cell_list.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include <cmath>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Lj_cell_list::Lj_cell_list(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
  }

  bool Lj_cell_list::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "cutoff"))
      {
        GET_OR_CHOOSE_A_REAL(cutoff, "", "")
        if (cutoff < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Force field cutoff have to non-negative.");
      }
      else if (string_cmp(t, "epsilon"))
      {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT(epsilon)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Epsilon have to be non-negative.");
      }
      else if (string_cmp(t, "sigma"))
      {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT(sigma)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Sigma have to be non-negative.");
      }
      else if (string_cmp(t, "set_neighborlist") || string_cmp(t, "neighborlist"))
      {
        FIND_OBJECT_BY_NAME(neighborlist, it)
        neighborlist = object_container->neighborlist[it->second.index];
      }
      else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Lj_cell_list::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(neighborlist)
  }

  void Lj_cell_list::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    auto cutoff_sq = cutoff * cutoff;

    const auto &pos = atom_data->owned.position;
    const auto pos_size = pos.size();
    const auto &binlist = neighborlist->binlist;
    const auto &nb = neighborlist->neigh_bin;

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < pos.size(); ++i)
    {
      const auto &pos_i = atom_data->owned.position[i];
      const auto type_i = atom_data->owned.type[i];
      const auto mass_inv_i = atom_data->owned.mass_inv[type_i];

      const auto nb_i = neighborlist->neigh_bin_index(pos_i);

      for (unsigned nb_j = 0; nb_j < nb[nb_i].size(); ++nb_j)
      {
        const auto &nb_ij = nb[nb_i][nb_j];

        for (unsigned nb_k = 0; nb_k < binlist[nb_ij.x][nb_ij.y][nb_ij.z].size(); ++nb_k)
        {

          unsigned j = binlist[nb_ij.x][nb_ij.y][nb_ij.z][nb_k];

          bool is_ghost = j >= pos_size;
          Vector<Real_t> pos_j;
          Real_t type_j, mass_inv_j;
          if (is_ghost)
          {
            j -= pos_size;
            pos_j = atom_data->ghost.position[j];
            type_j = atom_data->ghost.type[j];
          }
          else
          {
            pos_j = atom_data->owned.position[j];
            type_j = atom_data->owned.type[j];
            if (!(i < j))
              continue; // XXX CELL_LIST additional condition XXX//
          }
          mass_inv_j = atom_data->owned.mass_inv[type_j];

          auto dr = pos_j - pos_i;
          auto dr_sq = dr * dr;
          if (dr_sq > cutoff_sq)
            continue;
          const auto eps_ij = epsilon[type_i][type_j];
          const auto sigma_ij = sigma[type_i][type_j];

          auto r_c_sq_inv = 1 / (cutoff * cutoff); // XXX revise this after cutoff_list_activated addition
          auto rho_c_sq_inv = sigma_ij * sigma_ij * r_c_sq_inv;
          auto rho_c_6_inv = rho_c_sq_inv * rho_c_sq_inv * rho_c_sq_inv;
          auto rho_c_12_inv = rho_c_6_inv * rho_c_6_inv;
          auto dr_sq_inv = 1 / dr_sq;
          auto rho_sq_inv = sigma_ij * sigma_ij * dr_sq_inv;
          auto rho_6_inv = rho_sq_inv * rho_sq_inv * rho_sq_inv;
          auto rho_12_inv = rho_6_inv * rho_6_inv;

          auto force = 4 * eps_ij * (-12 * rho_12_inv * dr_sq_inv + 6 * rho_6_inv * dr_sq_inv + +12 * rho_c_12_inv * r_c_sq_inv - 6 * rho_c_6_inv * r_c_sq_inv) * dr;

          atom_data->owned.acceleration[i] += force * mass_inv_i;
          if (!is_ghost)
          {
#ifdef CAVIAR_WITH_OPENMP
#pragma omp atomic
            atom_data->owned.acceleration[j].x -= force.x * mass_inv_j;
#pragma omp atomic
            atom_data->owned.acceleration[j].y -= force.y * mass_inv_j;
#pragma omp atomic
            atom_data->owned.acceleration[j].z -= force.z * mass_inv_j;
#else
            atom_data->owned.acceleration[j] -= force * mass_inv_j;
#endif
          }
        }
      }
    }
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
