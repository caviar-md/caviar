
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

#include "caviar/objects/force_field/electrostatic_ewald_r.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/macro_constants.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Electrostatic_ewald_r::Electrostatic_ewald_r(CAVIAR *fptr) : Force_field{fptr},
                                                               k_electrostatic{1.0}
  {
    FC_OBJECT_INITIALIZE_INFO
    alpha = 1.0;
  }

  bool Electrostatic_ewald_r::read(caviar::interpreter::Parser *parser)
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
      else if (string_cmp(t, "alpha"))
      {
        GET_OR_CHOOSE_A_REAL(alpha, "", "")
        if (alpha < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "alpha have to non-negative.");
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
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Electrostatic_ewald_r::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(neighborlist)
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

  void Electrostatic_ewald_r::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    // XXX Working scheme (neighlist) with both 'cell_list' and 'verlet_list'
    ///*
    const auto &pos = atom_data->atom_struct_owned.position;
    const unsigned pos_size = pos.size();
    const auto alpha_sq = alpha * alpha;
    double virialLocal = 0;
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
        double lambda_coef = 1.0;

        if (lambda_is_set)
        {
          if (mol_index_i == mol_index_j)
          {
            if (lambda[type_i][type_j] == 0.0)
              continue;
            lambda_coef *= lambda[type_i][type_j];
          }
        }

        const auto r_ij = pos_j - pos_i;

        if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0)
          continue;

        const auto charge_j = atom_data->atom_type_params.charge[type_j];


        const auto rijml = r_ij;
        const auto rijml_sq = rijml * rijml;
        const auto rijml_norm = std::sqrt(rijml_sq);
        const auto erfc_arg = alpha * rijml_norm;

        // if (erfc_arg > erfc_cutoff) continue;

        // const auto sum_r = (2 * alpha * FC_PIS_INV * std::exp(-alpha_sq * rijml_sq) + std::erfc(erfc_arg) / rijml_norm) * (rijml / rijml_sq);

        // const auto force = - k_electrostatic * charge_i * charge_j * sum_r;

        const auto sum_r_x = (2 * alpha * FC_PIS_INV * std::exp(-alpha_sq * rijml_sq) + std::erfc(erfc_arg) / rijml_norm) * (1.0 / rijml_sq);

        auto forceCoef = -lambda_coef * k_electrostatic * charge_i * charge_j * sum_r_x;

        ////

        if (id_i < id_j)
        {
          virialLocal += -forceCoef * rijml_sq;
        }
        auto force = forceCoef * rijml;

        ////

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
    }
    //*/

    // XXX Working scheme (binlist) of  'cell_list'
    /*
      const auto &pos = atom_data -> atom_struct_owned.position;
      const unsigned pos_size = pos.size();
      const auto alpha_sq = alpha*alpha;

      const auto &binlist = neighborlist -> binlist;
      const auto &nb = neighborlist -> neigh_bin;


      for (unsigned i = 0; i < pos_size; ++i) {
        const auto pos_i = atom_data->atom_struct_owned.position [i];
        const auto type_i = atom_data -> atom_struct_owned.type [ i ];
        const auto charge_i = atom_data -> atom_type_params.charge [ type_i ];
        const auto mass_inv_i = atom_data -> atom_type_params.mass_inv [ type_i ];

        const auto nb_i = neighborlist -> neigh_bin_index (pos_i);
        for (unsigned nb_j = 0; nb_j < nb[nb_i].size(); ++nb_j) {
        const auto &nb_ij = nb[nb_i][nb_j];

        for (unsigned bl_i = 0; bl_i < binlist [nb_ij.x] [nb_ij.y] [nb_ij.z].size(); ++bl_i) {

          unsigned int j = binlist[nb_ij.x] [nb_ij.y] [nb_ij.z][bl_i];

          if (i>j) continue; // Works as multiplying a '0.5' to the sum. ...
          // ... We have to do this because  we use Newton's third law down there.
          //

          bool is_ghost = j >= pos_size;
          Vector<Real_t> pos_j;
          Real_t type_j;
          if (is_ghost) {
            j -= pos_size;
            pos_j = atom_data->atom_struct_ghost.position [j];
            type_j = atom_data->atom_struct_ghost.type [j];
          } else {
            pos_j = atom_data->atom_struct_owned.position [j];
            type_j = atom_data->atom_struct_owned.type [j];

          }

          const auto charge_j = atom_data -> atom_type_params.charge [ type_j ];
          const auto mass_inv_j = atom_data -> atom_type_params.mass_inv [ type_j ];
          const auto r_ij = pos_i - pos_j;

          if (r_ij.x==0 && r_ij.y==0 && r_ij.z==0) continue;

          const auto rijml = r_ij;
          const auto rijml_sq = rijml*rijml;
          const auto rijml_norm = std::sqrt(rijml_sq);
          const auto erfc_arg = alpha*rijml_norm;

          //if (erfc_arg > erfc_cutoff) continue;

          const auto sum_r = (2*alpha*FC_PIS_INV*std::exp(-alpha_sq*rijml_sq)
                           + std::erfc(erfc_arg) / rijml_norm )*(rijml/rijml_sq);

          const auto force =  k_electrostatic * charge_i * charge_j *sum_r;

          atom_data -> atom_struct_owned.acceleration[i] += force * mass_inv_i;
          if (!is_ghost)
            atom_data -> atom_struct_owned.acceleration[j] -= force * mass_inv_j;

        }
        }
      }
    */

    // XXX Working Scheme using field functions. Only 'binlist' field works.
    /*
     const auto &pos = atom_data -> atom_struct_owned.position;
     const unsigned pos_size = pos.size();

     for (unsigned i = 0; i < pos_size; ++i) {
       const auto pos_i = atom_data->atom_struct_owned.position [i];
       const auto type_i = atom_data -> atom_struct_owned.type [ i ];
       const auto charge_i = atom_data -> atom_type_params.charge [ type_i ];
       const auto mass_inv_i = atom_data -> atom_type_params.mass_inv [ type_i ];

       const auto force = charge_i * field (pos_i); // Working (binlist)
   //    const auto force = charge_i * field (i); // XXX (neighlist) won't work
       atom_data -> atom_struct_owned.acceleration[i] += force * mass_inv_i;

     }

    */
    atom_data->virialForce += virialLocal;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
