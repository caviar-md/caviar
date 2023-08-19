
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

#include "caviar/objects/force_field/gravity.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Gravity::Gravity(CAVIAR *fptr) : Force_field{fptr}, k_gravity{1.0}, external_field{Vector<double>{0, 0, 0}}
  {
    FC_OBJECT_INITIALIZE_INFO
  }

  bool Gravity::read(caviar::interpreter::Parser *parser)
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
      else if (string_cmp(t, "external_field"))
      {
        GET_OR_CHOOSE_A_REAL_3D_VECTOR(external_field, "", "");
        // if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Epsilon have to be non-negative.");
      }
      else if (string_cmp(t, "k_gravity"))
      {
        GET_OR_CHOOSE_A_REAL(k_gravity, "", "")
        if (k_gravity < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "k_gravity has to be non-negative.");
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

  void Gravity::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    my_mpi_rank = atom_data->get_mpi_rank();
  }

  void Gravity::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    const auto &pos = atom_data->atom_struct_owned.position;
    {
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
      for (unsigned int i = 0; i < pos.size(); ++i)
      {
        #ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
        const auto type_i = atom_data->atom_struct_owned.type[i];
        const auto mass_inv_i = atom_data->atom_type_params.mass_inv[type_i];

        for (unsigned int j = i + 1; j < pos.size(); ++j)
        {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
          const auto type_j = atom_data->atom_struct_owned.type[j];
          const auto mass_inv_j = atom_data->atom_type_params.mass_inv[type_j];
          const auto dr = pos[j] - pos[i];
          const auto dr_sq = dr * dr;
          const auto dr_norm = std::sqrt(dr_sq);
          const auto force = k_gravity * mass_inv_i * mass_inv_j * dr / (dr_sq * dr_norm);

          atom_data->atom_struct_owned.acceleration[i] += force * mass_inv_i;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[j].x -= force.x * mass_inv_j;
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[j].y -= force.y * mass_inv_j;
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[j].z -= force.z * mass_inv_j;
#else
          atom_data->atom_struct_owned.acceleration[j] -= force * mass_inv_j;
#endif
        }
      }
    }
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
