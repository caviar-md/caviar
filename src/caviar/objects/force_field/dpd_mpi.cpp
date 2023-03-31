
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

#include "caviar/objects/force_field/dpd_mpi.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"
#include <cmath>
#include <map>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Dpd_mpi::Dpd_mpi(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    kBoltzman = 1.0;
    temperature = 1.0;
    domain = nullptr;
  }

  bool Dpd_mpi::read(caviar::interpreter::Parser *parser)
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
      else if (string_cmp(t, "set_domain") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Dpd_mpi::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(neighborlist)
  }

  void Dpd_mpi::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

#ifdef CAVIAR_WITH_MPI
    auto neighborlist_domains = domain->neighborlist_domains;
    std::vector<int> g_num_recv, g_num_send;
    std::vector<std::vector<Vector<Real_t>>> g_send_accel, g_recv_accel;
    std::vector<std::vector<GlobalID_t>> g_send_id, g_recv_id;
    int nd = domain->neighborlist_domains.size();
    g_send_accel.resize(nd);
    g_recv_accel.resize(nd);
    g_send_id.resize(nd);
    g_recv_id.resize(nd);
    g_num_recv.resize(nd, 0);
    g_num_send.resize(nd, 0);

    std::map<int, int> rank_to_index;
    for (unsigned i = 0; i < domain->neighborlist_domains.size(); ++i)
    {                                             // this doesn't need to be reconstructed in every step
      rank_to_index[neighborlist_domains[i]] = i; // only after send_owned() happened it is needed to be cleared.
    }
#else
    std::vector<Vector<Real_t>> g_send_accel;
    std::vector<GlobalID_t> g_send_id;
#endif

    std::map<GlobalID_t, GlobalID_t> id_to_index;
    for (unsigned i = 0; i < atom_data->owned.id.size(); ++i)
    {                                          // this doesn't need to be reconstructed in every step
      id_to_index[atom_data->owned.id[i]] = i; // only after send_owned() happened it is needed to be cleared.
    }

    auto cutoff_sq = cutoff * cutoff;
    const auto dt_sq_inv = 1.0 / std::sqrt(dt);
    const auto &nlist = neighborlist->neighlist;
    for (unsigned int i = 0; i < nlist.size(); ++i)
    {
      const auto &pos_i = atom_data->owned.position[i];
      const auto &vel_i = atom_data->owned.velocity[i];
      const auto type_i = atom_data->owned.type[i];
      const auto mass_inv_i = atom_data->owned.mass_inv[type_i];
      const auto id_i = atom_data->owned.id[i];
      for (auto j : nlist[i])
      {
        bool is_ghost = j >= nlist.size();
        Vector<Real_t> pos_j, vel_j;
        Real_t type_j, mass_inv_j;
        GlobalID_t id_j = 0;
        if (is_ghost)
        {
          j -= nlist.size();
          id_j = atom_data->ghost.id[j];
          if (id_i > id_j)
          {
            continue; // there's another continue below which can make problems if we do a "++g_num_recv;" here
          }
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
        auto force = (force_conserv + force_dissip + force_rand) * dr_norm;

        atom_data->owned.acceleration[i] += force * mass_inv_i;
        if (!is_ghost)
          atom_data->owned.acceleration[j] -= force * mass_inv_j;
        else
        {
#ifdef CAVIAR_WITH_MPI
          int rti = rank_to_index[atom_data->ghost_rank[j]];
          g_send_accel[rti].push_back(-force * mass_inv_j);
          g_send_id[rti].push_back(id_j);
          g_num_send[rti]++;
#else
          g_send_accel.push_back(-force * mass_inv_j);
          g_send_id.push_back(id_j);
#endif
        }
      }
    }

#ifdef CAVIAR_WITH_MPI
    for (auto i = 1; i < nd; ++i)
    {
      MPI_Send(&g_num_send[i], 1, MPI_INT, neighborlist_domains[i], 0, mpi_comm);
      MPI_Recv(&g_num_recv[i], 1, MPI_INT, neighborlist_domains[i], 0, mpi_comm, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(mpi_comm);

    for (auto i = 1; i < nd; ++i)
    {
      int send_size = g_num_send[i];
      if (send_size != 0)
      {
        MPI_Send(&g_send_accel[i][0], 3 * send_size, MPI_DOUBLE, neighborlist_domains[i], 0, mpi_comm);
        MPI_Send(&g_send_id[i][0], send_size, MPI_INT, neighborlist_domains[i], 1, mpi_comm);
      }
    }

    for (auto i = 1; i < nd; ++i)
    {
      int recv_size = g_num_recv[i];
      if (recv_size != 0)
      {
        g_recv_id[i].resize(recv_size);
        g_recv_accel[i].resize(recv_size);
        MPI_Recv(&g_recv_accel[i][0], 3 * recv_size, MPI_DOUBLE, neighborlist_domains[i], 0, mpi_comm, MPI_STATUS_IGNORE);
        MPI_Recv(&g_recv_id[i][0], recv_size, MPI_INT, neighborlist_domains[i], 1, mpi_comm, MPI_STATUS_IGNORE);
      }
    }

    for (auto i = 1; i < nd; ++i)
    {
      for (auto j = 0; j < g_num_recv[i]; ++j)
      {
        int k = id_to_index.at(g_recv_id[i][j]);
        atom_data->owned.acceleration[k] += g_recv_accel[i][j];
      }
    }

    for (unsigned j = 0; j < g_send_id[0].size(); ++j)
    {
      int k = id_to_index.at(g_send_id[0][j]);
      atom_data->owned.acceleration[k] += g_send_accel[0][j]; // note that we didn't send this part.
    }
#else
    for (unsigned int j = 0; j < g_send_id.size(); ++j)
    {
      int k = id_to_index.at(g_send_id[j]);
      atom_data->owned.acceleration[k] += g_send_accel[j];
    }
#endif
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
