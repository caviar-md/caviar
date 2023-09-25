
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Implementation of Fene_bond is done by Nasrin Eyvazi
//
// This file is part of the CAVIAR package.
//
//
// The CAVIAR package is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the CAVIAR distribution.
//
//========================================================================

#include "caviar/objects/force_field/fene_bond.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Fene_bond::Fene_bond(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
  }

  bool Fene_bond::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "elastic_coef"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(elastic_coef)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Elastic coef. have to be non-negative.");
      }
      else if (string_cmp(t, "dissip_coef"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(dissip_coef)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Dissipation coef. have to be non-negative.");
      }
      else if (string_cmp(t, "epsilon_coef"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(epsilon_coef)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Epsilon coef. have to be non-negative.");
      }
      else if (string_cmp(t, "sigma_coef"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(sigma_coef)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Sigma coef. have to be non-negative.");
      }
      else if (string_cmp(t, "R_0"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(R)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "R_0 have to be non-negative.");
      }
      else if (string_cmp(t, "repulsive_cutoff"))
      {
        GET_OR_CHOOSE_A_REAL(repulsive_cutoff, "", "")
        if (repulsive_cutoff < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Repulsive cutoff have to non-negative.");
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
      else if (string_cmp(t, "set_domain") || string_cmp(t, "domain"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Fene_bond::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(domain)
    my_mpi_rank = atom_data->get_mpi_rank();
  }

  void Fene_bond::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    auto &pos = atom_data->atom_struct_owned.position;
    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &type = atom_data->atom_struct_owned.type;
    auto &mass_inv = atom_data->atom_type_params.mass_inv;
    double virialLocal = 0;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < atom_data->molecule_struct_owned.size(); i++)
    {
      auto &atomic_bond_vector = atom_data->molecule_struct_owned[i].atomic_bond_vector;

      for (unsigned int j = 0; j < atomic_bond_vector.size(); j++)
      {
        int id_1 = atomic_bond_vector[j].id_1;
        int id_2 = atomic_bond_vector[j].id_2;
        int k1 = atom_data->atom_id_to_index[id_1];
        int k2 = atom_data->atom_id_to_index[id_2];
        int btype = atomic_bond_vector[j].type;
        //        double d = atomic_bond_vector[i][j].length;
        //	double c = repulsive_cutoff * repulsive_cutoff;

#if defined(CAVIAR_WITH_MPI)
        const auto dr = pos[k2] - pos[k1];

#else
        const auto dr = domain->periodic_distance(pos[k2] - pos[k1]);
#endif
        const auto dv = vel[k2] - vel[k1];

        const auto dr_sq = dr * dr;
        const auto dr_norm = std::sqrt(dr_sq);
        const auto dr_vec = dr / dr_norm;
        auto cutoff_sq = repulsive_cutoff * repulsive_cutoff;
        auto c = cutoff_sq;
        double res;
        if (dr_sq <= c)
          res = 1.0;
        if (dr_sq > c)
          res = 0.0;
        const auto repulsive_term = 24 * (epsilon_coef[btype] / dr_norm) * (2 * pow(sigma_coef[btype] / dr_norm, 12) - pow(sigma_coef[btype] / dr_norm, 6)) * dr_vec;
        const auto force = -elastic_coef[btype] * (dr_norm / (1 - (dr_norm / R[btype]) * (dr_norm / R[btype]))) * dr_vec + repulsive_term * res - (dissip_coef[btype] * dv);

#ifdef CAVIAR_WITH_OPENMP
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k1].x -= force.x * mass_inv[type[k1]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k1].y -= force.y * mass_inv[type[k1]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k1].z -= force.z * mass_inv[type[k1]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k2].x += force.x * mass_inv[type[k2]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k2].y += force.y * mass_inv[type[k2]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k2].z += force.z * mass_inv[type[k2]];
#else
        atom_data->atom_struct_owned.acceleration[k1] -= force * mass_inv[type[k1]];
        atom_data->atom_struct_owned.acceleration[k2] += force * mass_inv[type[k2]];
#endif
      }
    }
    atom_data->virialForce += virialLocal;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
