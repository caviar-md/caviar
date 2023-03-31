
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

#include "caviar/objects/force_field/dpd.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include <cmath>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Dpd::Dpd(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO

    kBoltzman = 1.0;
    temperature = 1.0;
    dt = -1.0;
  }

  bool Dpd::read(caviar::interpreter::Parser *parser)
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
      else if (string_cmp(t, "temperature"))
      {
        GET_OR_CHOOSE_A_REAL(temperature, "", "")
        if (temperature < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "DPD temperature have to be non-negative.");
      }
      else if (string_cmp(t, "dt"))
      {
        GET_OR_CHOOSE_A_REAL(dt, "", "")
        if (dt < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "DPD dt have to be non-negative.");
      }
      else if (string_cmp(t, "kb"))
      {
        GET_OR_CHOOSE_A_REAL(kBoltzman, "", "")
        if (kBoltzman < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "DPD kBoltzman have to be non-negative.");
      }
      else if (string_cmp(t, "seed") || string_cmp(t, "random_seed"))
      {
        GET_OR_CHOOSE_A_REAL(rnd_seed, "", "")
        if (rnd_seed < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "random_seed have to non-negative.");
        rnd_generator.seed(rnd_seed);
      }
      else if (string_cmp(t, "conserv_coef"))
      {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT(conserv_coef)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "conserv_coef have to be non-negative.");
      }
      else if (string_cmp(t, "dissip_coef"))
      {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT(dissip_coef)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "dissip_coef have to be non-negative.");
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

  void Dpd::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(neighborlist)
  }

  void Dpd::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    auto cutoff_sq = cutoff * cutoff;
    auto dt_sq_inv = 1.0 / std::sqrt(dt);
    const auto &nlist = neighborlist->neighlist;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < nlist.size(); ++i)
    {
      const auto &pos_i = atom_data->owned.position[i];
      const auto &vel_i = atom_data->owned.velocity[i];
      const auto type_i = atom_data->owned.type[i];
      const auto mass_inv_i = atom_data->owned.mass_inv[type_i];
      for (auto j : nlist[i])
      {
        bool is_ghost = j >= nlist.size();
        Vector<Real_t> pos_j, vel_j;
        Real_t type_j, mass_inv_j;
        if (is_ghost)
        {
          j -= nlist.size();
          pos_j = atom_data->ghost.position[j];
          vel_j = atom_data->ghost.velocity[j];
          type_j = atom_data->ghost.type[j];
        }
        else
        {
          pos_j = atom_data->owned.position[j];
          vel_j = atom_data->owned.velocity[j];
          type_j = atom_data->owned.type[j];
        }
        mass_inv_j = atom_data->owned.mass_inv[type_j];
        auto dr = pos_j - pos_i;
        auto dv = vel_j - vel_i;
        auto r_sq = dr * dr;
        auto r_sqrt = std::sqrt(r_sq);
        if (r_sq > cutoff_sq)
          continue;
        auto dr_norm = dr / r_sqrt;
        const auto conserv_coef_ij = conserv_coef[type_i][type_j];
        const auto dissip_coef_ij = dissip_coef[type_i][type_j];
        auto sigma = std::sqrt(2.0 * kBoltzman * temperature * dissip_coef_ij);
        auto w_r = 1.0 - (r_sqrt / cutoff);
        auto alpha = rnd_ndist(rnd_generator); //
        auto force_conserv = conserv_coef_ij * w_r;
        auto force_dissip = -dissip_coef_ij * w_r * w_r * (dr_norm * dv);
        auto force_rand = sigma * w_r * alpha * dt_sq_inv;
        auto force = -(force_conserv + force_dissip + 0.0 * force_rand) * dr_norm;
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

} // force_field

CAVIAR_NAMESPACE_CLOSE
