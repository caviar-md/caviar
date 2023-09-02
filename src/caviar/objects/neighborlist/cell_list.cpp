
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

#include "caviar/objects/neighborlist/cell_list.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/domain.h"

#include <cmath>

CAVIAR_NAMESPACE_OPEN

namespace neighborlist
{
  inline int int_floor(double x)
  {
    return (int)(x + 100000) - 100000;
  }

  void remove_duplicates(std::vector<Vector<int>> &v)
  {
    auto end = v.end();
    for (auto it = v.begin(); it != end; ++it)
    {
      end = std::remove(it + 1, end, *it);
    }

    v.erase(end, v.end());
  }

  Cell_list::Cell_list(CAVIAR *fptr) : Neighborlist{fptr}, domain{nullptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    make_neighlist = true;
  }

  bool Cell_list::read(caviar::interpreter::Parser *parser)
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
        if (cutoff <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "cutoff have to non-negative.");
      }
      else if (string_cmp(t, "cutoff_extra_coef"))
      {
        GET_OR_CHOOSE_A_REAL(cutoff_extra_coef, "", "")
        if (cutoff_extra_coef < 1.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "cutoff_extra_coef must be larger than 1.0");
      }
      else if (string_cmp(t, "make_neighlist"))
      {
        make_neighlist = true;
      }
      else if (string_cmp(t, "set_domain") || string_cmp(t, "domain"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
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

  void Cell_list::init()
  {
    if (!initialize)
      return;
    initialize = false;

    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(domain)

    if (cutoff <= 0.0)
      error->all(FC_FILE_LINE_FUNC, "cutoff have to non-negative.");
    if (cutoff_extra_coef < 1.0)
      error->all(FC_FILE_LINE_FUNC, "cutoff_extra_coef must be equal or larger than 1.0");

    cutoff_extra = cutoff * cutoff_extra_coef;
    threshold_distance_sq = 0.25 * (cutoff_extra - cutoff) * (cutoff_extra - cutoff);
    cutoff_inv = 1.0 / cutoff_extra;

    if (!make_neighlist)
    {
      output->warning("'make_neighlist' is false. The force_fields using verlet_list won't function.");
    }

    my_mpi_rank = atom_data->get_mpi_rank();

    auto dd = domain->upper_local - domain->lower_local; // NO DOMAIN DECOMPOSITION


    no_bins.x = std::ceil(dd.x * cutoff_inv);
    no_bins.y = std::ceil(dd.y * cutoff_inv);
    no_bins.z = std::ceil(dd.z * cutoff_inv);

    binlist.resize(no_bins.x);

    for (unsigned i = 0; i < binlist.size(); ++i)
      binlist[i].resize(no_bins.y);

    for (unsigned i = 0; i < binlist.size(); ++i)
      for (unsigned j = 0; j < binlist[i].size(); ++j)
        binlist[i][j].resize(no_bins.z);

    make_neigh_bin();


    // std::cout << "   " << cutoff << std::endl;
    // std::cout << "cutoff : " << cutoff << std::endl;
    // std::cout << "cutoff_extra : " << cutoff_extra << std::endl;
    // std::cout << "threshold_distance_sq : " << threshold_distance_sq << std::endl;
    // std::cout << "cutoff_inv : " << cutoff_inv << std::endl;

    // // std::cout << "domain->upper_local : " << domain->upper_local << std::endl;
    // // std::cout << "domain->lower_local : " << domain->lower_local << std::endl;

    // // std::cout << "dd : " << dd << std::endl;
    // // std::cout << "no_bins : " << no_bins << std::endl;
    // // std::cout << "binlist_index(0,0,0)  : " << binlist_index(Vector<double>(0, 0, 0)) << std::endl;
    // // std::cout << "binlist_index(-120,-12,-12)  : " << binlist_index(Vector<double>(-120, -12, -12)) << std::endl;
    // // std::cout << "binlist_index(120,12,12)  : " << binlist_index(Vector<double>(120, 12, 12)) << std::endl;
    // for (int i = 0; i < no_bins.x; ++i)
    // {
    //   for (int j = 0; j < no_bins.y; ++j)
    //   {
    //     for (int k = 0; k < no_bins.z; ++k)
    //     {
    //       int m = i + no_bins.x * j + no_bins.x * no_bins.y * k;
    //       std::cout << "======================" << std::endl;
 
    //       std::cout << m << " : [" << i << "," << j << "," << k << "] : ( " << (i *cutoff_extra) + domain->lower_local.x << " , " << (j / cutoff_inv) + domain->lower_local.y << " , " << (k / cutoff_inv) + domain->lower_local.z << " )" << std::endl;

    //       // for (int n = 0; n < neigh_bin[m].size(); ++n)
    //       // {
    //       //   std::cout << " , [" << neigh_bin[m][n].x << "," << neigh_bin[m][n].y << "," << neigh_bin[m][n].z << "] ";
    //       // }
    //       // std::cout << std::endl;
    //     }
    //   }
    // }
  }

  Vector<int> Cell_list::binlist_index(const Vector<double> &pos)
  {
    Vector<int> ind;
    ind.x = std::floor((pos.x - domain->lower_local.x) * cutoff_inv); // int_floor
    ind.y = std::floor((pos.y - domain->lower_local.y) * cutoff_inv);
    ind.z = std::floor((pos.z - domain->lower_local.z) * cutoff_inv);
    //================================================
    // Handling Ghosts in periodic boundary condition
    //================================================
    if (ind.x < 0)
      ind.x = 0;
    if (ind.y < 0)
      ind.y = 0;
    if (ind.z < 0)
      ind.z = 0;
    if (ind.x >= no_bins.x)
      ind.x = no_bins.x - 1;
    if (ind.y >= no_bins.y)
      ind.y = no_bins.y - 1;
    if (ind.z >= no_bins.z)
      ind.z = no_bins.z - 1;
    return ind;
  }

  int Cell_list::neigh_bin_index(const Vector<double> &pos)
  {
    auto ind = binlist_index(pos);
    return ind.x + no_bins.x * ind.y + no_bins.x * no_bins.y * ind.z;
  }

  void Cell_list::build_binlist()
  {
    for (unsigned int i = 0; i < binlist.size(); ++i)
    {
      for (unsigned int j = 0; j < binlist[i].size(); ++j)
      {
        for (unsigned int k = 0; k < binlist[i][j].size(); ++k)
        {
          binlist[i][j][k].clear();
        }
      }
    }

    const auto &pos = atom_data->atom_struct_owned.position;
    const auto &pos_ghost = atom_data->atom_struct_ghost.position;
    const auto pos_size = pos.size();
    const auto ghost_size = pos_ghost.size();

    for (unsigned int i = 0; i < pos_size; ++i)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      auto ind = binlist_index(pos[i]);
      binlist[ind.x][ind.y][ind.z].emplace_back(i);
    }

    for (unsigned int i = 0; i < ghost_size; ++i)
    {
      auto ind = binlist_index(pos_ghost[i]);
      binlist[ind.x][ind.y][ind.z].emplace_back(i + pos_size);
    }
  }

  void Cell_list::build_verlet_list_from_bins()
  {
    if (all_atom_test)
    {
      all_atom_test_function(0);
      return;
    }

    const auto &pos = atom_data->atom_struct_owned.position;
    const auto &pos_ghost = atom_data->atom_struct_ghost.position;

    const int pos_size = pos.size();
    const auto cutoff_extra_sq = cutoff_extra * cutoff_extra;
    const auto &nb = neigh_bin;
    neighlist.resize(pos_size);

    for (int i = 0; i < pos_size; ++i)
    {
      neighlist[i].clear();

#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif

      auto nb_i = neigh_bin_index(pos[i]);
      for (unsigned nb_j = 0; nb_j < nb[nb_i].size(); ++nb_j)
      {
        const auto &nb_ij = nb[nb_i][nb_j];
        for (unsigned j = 0; j < binlist[nb_ij.x][nb_ij.y][nb_ij.z].size(); ++j)
        {
          auto dr = pos[i];
          auto ind_j = binlist[nb_ij.x][nb_ij.y][nb_ij.z][j];
          bool ghost = (ind_j >= pos_size);
          if (ghost)
          {
            dr -= pos_ghost[ind_j - pos_size];
          }
          else
          {
            if (ind_j <= i)
              continue;

            dr -= pos[ind_j];
          }

          if (dr * dr < cutoff_extra_sq)
            neighlist[i].emplace_back(ind_j);
        }
      }
    }

    pos_old = pos;
    ghost_pos_old = pos_ghost;
#ifdef CAVIAR_WITH_MPI
    mpi_rank_old = atom_data->atom_struct_owned.mpi_rank;
#endif
  }

  void Cell_list::build_neighlist()
  {

    build_binlist();
    if (make_neighlist)
    {
      if (rebuild_neighlist())
      {
        build_verlet_list_from_bins();
      }
    }
  }

  // this function will be called only once at the start or after any change in
  // the Cell_list.
  // Here we make the 'neigh_bin'. It contains the list of neighbors (of cell type)
  // any cell
  // could have, including itself.
  void Cell_list::make_neigh_bin()
  {
    neigh_bin.clear();
    neigh_bin.resize(no_bins.x * no_bins.y * no_bins.z);
    auto bc = domain->boundary_condition;

    Vector<int> ind{0, 0, 0};

    for (ind.x = 0; ind.x < no_bins.x; ++ind.x)
    {
      for (ind.y = 0; ind.y < no_bins.y; ++ind.y)
      {
        for (ind.z = 0; ind.z < no_bins.z; ++ind.z)
        {
          int ind_x_min = ind.x - 1;
          int ind_x_max = ind.x + 1;
          int ind_y_min = ind.y - 1;
          int ind_y_max = ind.y + 1;
          int ind_z_min = ind.z - 1;
          int ind_z_max = ind.z + 1;

#if defined(CAVIAR_WITH_MPI) && !defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
          if (ind_x_min == -1)
            ind_x_min = 0;
          if (ind_x_max == no_bins.x)
            ind_x_max = no_bins.x - 1;

          if (ind_y_min == -1)
            ind_y_min = 0;
          if (ind_y_max == no_bins.y)
            ind_y_max = no_bins.y - 1;

          if (ind_z_min == -1)
            ind_z_min = 0;
          if (ind_z_max == no_bins.z)
            ind_z_max = no_bins.z - 1;
#else
          if (bc.x == 0)
          {
            if (ind_x_min == -1)
              ind_x_min = 0;
            if (ind_x_max == no_bins.x)
              ind_x_max = no_bins.x - 1;
          }

          if (bc.y == 0)
          {
            if (ind_y_min == -1)
              ind_y_min = 0;
            if (ind_y_max == no_bins.y)
              ind_y_max = no_bins.y - 1;
          }

          if (bc.z == 0)
          {
            if (ind_z_min == -1)
              ind_z_min = 0;
            if (ind_z_max == no_bins.z)
              ind_z_max = no_bins.z - 1;
          }

#endif

          int k = ind.x + no_bins.x * ind.y + no_bins.x * no_bins.y * ind.z;

          for (int ind_x_i = ind_x_min; ind_x_i <= ind_x_max; ++ind_x_i)
          {
            for (int ind_y_i = ind_y_min; ind_y_i <= ind_y_max; ++ind_y_i)
            {
              for (int ind_z_i = ind_z_min; ind_z_i <= ind_z_max; ++ind_z_i)
              {
#if defined(CAVIAR_WITH_MPI) && !defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
                caviar::Vector<int> nvec{ind_x_i, ind_y_i, ind_z_i};

                neigh_bin[k].emplace_back(nvec);

#else
                caviar::Vector<int> nvec{ind_x_i, ind_y_i, ind_z_i};

                if (nvec.x == no_bins.x)
                  nvec.x = 0;
                if (nvec.y == no_bins.y)
                  nvec.y = 0;
                if (nvec.z == no_bins.z)
                  nvec.z = 0;
                if (nvec.x == -1)
                  nvec.x = no_bins.x - 1;
                if (nvec.y == -1)
                  nvec.y = no_bins.y - 1;
                if (nvec.z == -1)
                  nvec.z = no_bins.z - 1;

                neigh_bin[k].emplace_back(nvec);
#endif
              }
            }
          }
        }
      }
    }

    for (unsigned int i = 0; i < neigh_bin.size(); ++i)
    {
      remove_duplicates(neigh_bin[i]);
    }
  }
} // neighborlist

CAVIAR_NAMESPACE_CLOSE
