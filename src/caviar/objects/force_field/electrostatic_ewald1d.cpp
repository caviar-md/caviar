
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

#include "caviar/objects/force_field/electrostatic_ewald1d.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Electrostatic_ewald1d::Electrostatic_ewald1d(CAVIAR *fptr) : Force_field{fptr}, k_electrostatic{1.0}
  {
    FC_OBJECT_INITIALIZE_INFO
    num_mirrors = 1;
    sigma = 1.0;
  }

  bool Electrostatic_ewald1d::read(caviar::interpreter::Parser *parser)
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
      else if (string_cmp(t, "num_mirrors"))
      {
        GET_OR_CHOOSE_A_INT(num_mirrors, "", "")
        if (num_mirrors < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "num_mirrors has to be non-negative.");
      }
      else if (string_cmp(t, "sigma"))
      {
        GET_OR_CHOOSE_A_REAL(sigma, "", "")
        if (sigma < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "sigma has to be non-negative.");
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

  void Electrostatic_ewald1d::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(domain)
    FC_NULLPTR_CHECK(neighborlist)
    my_mpi_rank = atom_data->get_mpi_rank();

    auto dd = domain->upper_local - domain->lower_local; // NO DOMAIN DECOMPOSITION
    auto bc = domain->boundary_condition;

    if (bc.x + bc.y + bc.z != 1)
      error->all(FC_FILE_LINE_FUNC, "expected one periodicity in the domain.");

    Vector<double> v{0, 0, 0};

    if (bc.x == 1)
      v.x = dd.x;
    if (bc.y == 1)
      v.y = dd.y;
    if (bc.z == 1)
      v.z = dd.z;

    lattice_vec.resize(2 * num_mirrors + 1);
    for (unsigned int i = 0; i < lattice_vec.size(); ++i)
    {
      int j = -num_mirrors + i;
      lattice_vec[i] = j * v;
      // std::cout <<lattice_vec[i] << "\n";
    }

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

  void Electrostatic_ewald1d::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    // XXX scheme (neighlist) with both 'cell_list' and 'verlet_list'
    ///*
    const auto &pos = atom_data->atom_struct_owned.position;
    const unsigned pos_size = pos.size();

    const auto lattice_vec_size = lattice_vec.size();
    const auto sigma_sq = sigma * sigma;
    double virialLocal = 0;

    bool get_pressure_process = atom_data->get_pressure_process();

    const auto &nlist = neighborlist->neighlist;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < pos_size; ++i)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      const auto pos_i = atom_data->atom_struct_owned.position[i];
      const auto type_i = atom_data->atom_struct_owned.type[i];
      const auto charge_i = atom_data->atom_type_params.charge[type_i];
      const auto mass_inv_i = atom_data->atom_type_params.mass_inv[type_i];
      int id_i = atom_data->atom_struct_owned.id[i];
      const auto mol_index_i = atom_data->atom_struct_owned.molecule_index[i];

      // short range part
      for (auto j : nlist[i])
      {
        bool is_ghost = j >= pos_size;
        Vector<Real_t> pos_j;
        Real_t type_j;
        int id_j;
        int mol_index_j;

        if (is_ghost)
        {
          j -= pos_size;
          pos_j = atom_data->atom_struct_ghost.position[j];
          type_j = atom_data->atom_struct_ghost.type[j];
          id_j = atom_data->atom_struct_ghost.id[j];
          mol_index_j = atom_data->atom_struct_ghost.molecule_index[j];
        }
        else
        {
          pos_j = atom_data->atom_struct_owned.position[j];
          type_j = atom_data->atom_struct_owned.type[j];
          id_j = atom_data->atom_struct_owned.id[j];
          mol_index_j = atom_data->atom_struct_owned.molecule_index[j];
        }

        const auto charge_j = atom_data->atom_type_params.charge[type_j];
        const auto r_ij = pos_j - pos_i;

        if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0)
          continue;

        const auto r_ij_sq = r_ij * r_ij;
        const auto d1 = 1.0 / std::sqrt(r_ij_sq);
        const auto d2 = 1.0 / std::sqrt(r_ij_sq + sigma_sq);

        // const auto sum_r = r_ij * ((d1 * d1 * d1) - (d2 * d2 * d2));

        // const auto force = - k_electrostatic * charge_i * charge_j * sum_r;

        const auto sum_r_x = ((d1 * d1 * d1) - (d2 * d2 * d2));

        auto forceCoef = -k_electrostatic * charge_i * charge_j * sum_r_x;

        if (lambda_is_set)
        {
          if (mol_index_i == mol_index_j)
          {
            forceCoef *= lambda[type_i][type_j];
          }
        }

        if (id_i < id_j)
        {
          virialLocal += -forceCoef * r_ij_sq;
        }
        auto force = forceCoef * r_ij;

        atom_data->atom_struct_owned.acceleration[i] += force * mass_inv_i;
        if (!is_ghost)
        {
          const auto mass_inv_j = atom_data->atom_type_params.mass_inv[type_j];

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

      ///*
      // long range part

      Vector<double> field{0, 0, 0};
      for (unsigned int j = 0; j < pos.size(); ++j)
      {
#ifdef CAVIAR_WITH_MPI
        if (atom_data->atom_struct_owned.mpi_rank[j] != my_mpi_rank)
          continue;
#endif
        const auto type_j = atom_data->atom_struct_owned.type[j];
        const auto charge_j = atom_data->atom_type_params.charge[type_j];

        Vector<double> sum{0, 0, 0};
        for (unsigned int k = 0; k < lattice_vec_size; ++k)
        {
          const auto dr = pos[i] - pos[j] + lattice_vec[k];
          const auto dr_sq = dr * dr;
          // if (dr_sq == 0.0) continue; //XXX removing self field?
          const auto d2 = 1.0 / std::sqrt(dr_sq + sigma_sq);
          sum += dr * (d2 * d2 * d2);
        }
        field += charge_j * sum;
      }

      auto force = k_electrostatic * charge_i * field;
      // if (lambda_is_set) force *= lambda[type_i][type_i];

      if (get_pressure_process)
      {
        atom_data->add_to_external_virial(force, i);
        // or ???
        // atom_data->add_to_external_virial(force, i, p_i);
      }

      atom_data->atom_struct_owned.acceleration[i] += force * mass_inv_i;

      //*/
    }

    //*/

    // XXX Scheme using field functions.
    /*
      const auto &pos = atom_data -> atom_struct_owned.position;
      const unsigned pos_size = pos.size();

      for (unsigned i = 0; i < pos_size; ++i) {
        const auto pos_i = atom_data->atom_struct_owned.position [i];
        const auto type_i = atom_data -> atom_struct_owned.type [ i ];
        const auto charge_i = atom_data -> atom_type_params.charge [ type_i ];
        const auto mass_inv_i = atom_data -> atom_type_params.mass_inv [ type_i ];

        //const auto force = charge_i * field (pos_i);  //
        const auto force = charge_i * field (i);  //
        atom_data -> atom_struct_owned.acceleration[i] += force * mass_inv_i;

      }
    */

    atom_data->virialForce += virialLocal;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
