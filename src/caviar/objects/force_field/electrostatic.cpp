
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

#include "caviar/objects/force_field/electrostatic.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Electrostatic::Electrostatic(CAVIAR *fptr) : Force_field{fptr}, k_electrostatic{1.0}, external_field{Vector<double>{0, 0, 0}}
  {
    FC_OBJECT_INITIALIZE_INFO
  }

  bool Electrostatic::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR
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
      else if (string_cmp(t, "k_electrostatic"))
      {
        GET_OR_CHOOSE_A_REAL(k_electrostatic, "", "")
        if (k_electrostatic < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "k_electrostatic has to be non-negative.");
      }
      else if (string_cmp(t, "lambda"))
      {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT_R(lambda, 1.0)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Sigma have to be non-negative.");
        lambda_is_set = true;
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

  void Electrostatic::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(domain)
    my_mpi_rank = atom_data->get_mpi_rank();

    unsigned atom_data_type_max = 0;
    for (auto &&t : atom_data->atom_struct_owned.type)
    {
      if (atom_data_type_max < t)
        atom_data_type_max = t;
    }

    if (lambda_is_set)
    {
      if (lambda.size() <= atom_data_type_max)
        lambda.resize(atom_data_type_max + 1);
      for (unsigned int i = 0; i < lambda.size(); ++i)
      {
        if (lambda[i].size() <= atom_data_type_max)
          lambda[i].resize(atom_data_type_max + 1, 1.0);
      }
    }
  }

  void Electrostatic::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    double virialLocal = 0;

    bool get_pressure_process = atom_data->get_pressure_process();

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
        const auto charge_i = atom_data->atom_type_params.charge[type_i];
        int id_i = atom_data->atom_struct_owned.id[i];
        const auto mol_index_i = atom_data->atom_struct_owned.molecule_index[i];

        for (unsigned int j = i + 1; j < pos.size(); ++j)
        {

          const auto type_j = atom_data->atom_struct_owned.type[j];

          const auto mass_inv_j = atom_data->atom_type_params.mass_inv[type_j];
          const auto charge_j = atom_data->atom_type_params.charge[type_j];
          int id_j = atom_data->atom_struct_owned.id[i];

#if defined(CAVIAR_WITH_MPI)
          const auto dr = pos[j] - pos[i];
#else
          const auto dr = domain->periodic_distance(pos[j] - pos[i]);
#endif

          const auto dr_sq = dr * dr;
          const auto dr_norm = std::sqrt(dr_sq);
          auto forceCoef = -k_electrostatic * charge_i * charge_j / (dr_sq * dr_norm);
          const auto mol_index_j = atom_data->atom_struct_owned.molecule_index[j];

          if (lambda_is_set)
          {
            if (mol_index_i == mol_index_j)
            {
              forceCoef *= lambda[type_i][type_j];
            }
          }
          // const auto force = k_electrostatic * charge_i * charge_j * dr / (dr_sq * dr_norm);

          if (id_i < id_j)
          {
            virialLocal += -forceCoef * dr_sq;
          }
          auto force = forceCoef * dr;

          // std::cout << "force: " << force << std::endl;
          atom_data->atom_struct_owned.acceleration[i] += force * mass_inv_i;
          // std::cout << "acci: " << atom_data->atom_struct_owned.acceleration[i] << std::endl;

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
          // std::cout << "accj: " << atom_data->atom_struct_owned.acceleration[j] << std::endl;
        }

        auto force = external_field * charge_i;
        // if (lambda_is_set)
        //   force *= lambda[type_i][type_i];

        if (get_pressure_process)
        {
          atom_data->add_to_external_virial(force, i);
          // or ???
          // atom_data->add_to_external_virial(force, i, p_i);
        }

        atom_data->atom_struct_owned.acceleration[i] += force * mass_inv_i;
      }
    }
    atom_data->virialForce += virialLocal;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
