
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

#include "caviar/objects/neighborlist.h"
#include "caviar/interpreter/error.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

inline int int_floor(double x)
  {
    return (int)(x + 100000) - 100000;
  }

static void remove_duplicates(std::vector<Vector<int>> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = std::remove(it + 1, end, *it);
  }

  v.erase(end, v.end());
}

static void remove_duplicates(std::vector<unsigned int> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = std::remove(it + 1, end, *it);
  }

  v.erase(end, v.end());
}

Neighborlist::Neighborlist(CAVIAR *fptr) : Pointers{fptr},
                                           atom_data{nullptr} {
                                               FC_OBJECT_INITIALIZE}

                                           Neighborlist::~Neighborlist()
{
}

void Neighborlist::verify_settings()
{
}

bool Neighborlist::build(bool update_neighborlist)
{
  if (build_cell_list)
  {
    update_cell_list();
  }

  if (!update_neighborlist)
  {
    update_neighborlist = update_is_needed();
  }

  
  if (update_neighborlist)
  {
    if (make_verlet_list_from_cell_list)
      update_verlet_list_from_cell_list();
    else
      update_verlet_list();
    return true;
  }
  return false;
}


bool Neighborlist::update_is_needed()
{
  if (!initialized)
  {
    init();
  }

  if (rebuild_test)
  {
    return true;
  }

  const auto &pos = atom_data->atom_struct_owned.position;
  const auto &ghost_pos = atom_data->atom_struct_ghost.position;

  unsigned int pos_size = pos.size();
  unsigned int ghost_size = ghost_pos.size();

  if (pos_size != pos_old.size() || ghost_size != ghost_pos_old.size())
  {
    return true;
  }
  else
  {
    for (unsigned int i = 0; i < pos_size; ++i)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != mpi_rank_old[i]) // make new verlet_list if any mpi_rank is changed,
        return true;

      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank) // then, ignore mpi_ranks out of the domain
        continue;
#endif

      auto disp = pos[i] - pos_old[i];
      if (disp * disp > threshold_distance_sq)
        return true;
    }

    for (unsigned int i = 0; i < ghost_size; ++i)
    {
      auto disp = ghost_pos[i] - ghost_pos_old[i];
      if (disp * disp > threshold_distance_sq)
        return true;
    }
  }

  return false;
}

void Neighborlist::calculate_cutoff_extra()
{
  // if (cutoff_extra_coef <= 0)
  // {
  //   cutoff_extra = 0;
  //   return;
  // }

  // const auto &vel = atom_data->atom_struct_owned.velocity;

  // double max_vel_sq = 0.0;
  // for (unsigned int i = 0; i < vel.size(); ++i)
  // {
  //   double vel_sq_temp = vel[i] * vel[i];
  //   if (max_vel_sq < vel_sq_temp)
  //     max_vel_sq = vel_sq_temp;
  // }
  // cutoff_extra = cutoff_extra_coef * dt * std::sqrt(max_vel_sq);
}

void Neighborlist::all_atom_test_function(int state)
{
  if (state == 0)
  {
    const auto &pos = atom_data->atom_struct_owned.position;
    const auto &pos_ghost = atom_data->atom_struct_ghost.position;
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
        neighlist[i].emplace_back(j);
      }
      for (unsigned int j = 0; j < pos_ghost.size(); ++j)
      {
        neighlist[i].emplace_back(j + pos_size);
      }
    }
  }
}

bool Neighborlist::read(caviar::interpreter::Parser *parser)
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
    else if (string_cmp(t, "make_verlet_list_from_cell_list"))
    {
      make_verlet_list_from_cell_list = true;
      build_cell_list = true;
    }
      else if (string_cmp(t, "build_cell_list"))
    {
      build_cell_list = true;
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
    else if (string_cmp(t, "dt"))
    {
      GET_OR_CHOOSE_A_REAL(dt, "", "")
      if (dt < 0.0)
        error->all(FC_FILE_LINE_FUNC_PARSE, "dt have to non-negative.");
    }
    else
      error->all(FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
  }
  return in_file;
}

//============================================
//               Verlet_list functions
//============================================

void Neighborlist::update_verlet_list()
{

  if (all_atom_test)
  {
    all_atom_test_function(0);
    return;
  }

  const auto &pos = atom_data->atom_struct_owned.position;
  const auto &pos_ghost = atom_data->atom_struct_ghost.position;
  const auto cutoff_extra_sq = cutoff_extra * cutoff_extra;
  const auto pos_size = pos.size();
  neighlist.resize(pos_size);
  // std::cout << cutoff_extra_sq << std::endl;
  for (unsigned int i = 0; i < pos_size; ++i)
  {
    neighlist[i].clear();

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
      if (dr * dr < cutoff_extra_sq)
      {
        neighlist[i].emplace_back(j);
      }
    }
    for (unsigned int j = 0; j < pos_ghost.size(); ++j)
    {
      auto dr = pos_ghost[j] - pos[i];
      if (dr * dr < cutoff_extra_sq)
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

//============================================
//               cell_list functions
//============================================

void Neighborlist::init()
{
  FC_NULLPTR_CHECK(atom_data)

  my_mpi_rank = atom_data->get_mpi_rank();

  //============================================
  //               verlet_list init
  //============================================



  if (cutoff <= 0.0)
    error->all(FC_FILE_LINE_FUNC, "cutoff have to non-negative.");

  if (cutoff_extra_coef < 1.0)
    error->all(FC_FILE_LINE_FUNC, "cutoff_extra_coef must be equal or larger than 1.0");

  cutoff_extra = cutoff * cutoff_extra_coef;
  threshold_distance_sq = 0.25 * (cutoff_extra - cutoff) * (cutoff_extra - cutoff);

  const auto &pos = atom_data->atom_struct_owned.position;
  pos_old.resize(pos.size());

  //============================================
  //               cell_list init
  //============================================
  if (!build_cell_list)
    return;

  FC_NULLPTR_CHECK(domain)

  if (cutoff <= 0.0)
    error->all(FC_FILE_LINE_FUNC, "cutoff have to non-negative.");
  if (cutoff_extra_coef < 1.0)
    error->all(FC_FILE_LINE_FUNC, "cutoff_extra_coef must be equal or larger than 1.0");

  cutoff_extra = cutoff * cutoff_extra_coef;
  threshold_distance_sq = 0.25 * (cutoff_extra - cutoff) * (cutoff_extra - cutoff);
  cutoff_inv = 1.0 / cutoff_extra;

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

  // binlist_linear.resize(no_bins.x *no_bins.y * no_bins.z);
  
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

Vector<int> Neighborlist::binlist_index(const Vector<double> &pos)
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

int Neighborlist::neigh_bin_index(const Vector<double> &pos)
{
  auto ind = binlist_index(pos);
  return ind.x + no_bins.x * ind.y + no_bins.x * no_bins.y * ind.z;
}


void Neighborlist::update_cell_list()
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


/*

void Neighborlist::update_cell_list()
{
  for (unsigned int i = 0; i < binlist_linear.size(); ++i)
  {

    binlist_linear[i].clear();         
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
    binlist_linear[ind.x + no_bins.x * ind.y + no_bins.x * no_bins.y * ind.z].emplace_back(i);
  }

  for (unsigned int i = 0; i < ghost_size; ++i)
  {
    auto ind = binlist_index(pos_ghost[i]);
    binlist_linear[ind.x + no_bins.x * ind.y + no_bins.x * no_bins.y * ind.z].emplace_back(i + pos_size);
  }
}
*/


//==============================================================
// The result of this function (update_verlet_list_from_cell_lis) 
// is the same as (update_verlet_list)
// except it is not sorted. This makes different MD results 
// due to limited precision of floating point arithmetics
//==============================================================
void Neighborlist::update_verlet_list_from_cell_list()
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
          neighlist[i].emplace_back(ind_j); // not sorted
      }

    }
  }

  pos_old = pos;
  ghost_pos_old = pos_ghost;
#ifdef CAVIAR_WITH_MPI
  mpi_rank_old = atom_data->atom_struct_owned.mpi_rank;
#endif
}


// this function will be called only once at the start or after any change in
// the Cell_list.
// Here we make the 'neigh_bin'. It contains the list of neighbors (of cell type)
// any cell
// could have, including itself.
void Neighborlist::make_neigh_bin()
{
  neigh_bin.clear();
  neigh_bin.resize(no_bins.x * no_bins.y * no_bins.z);

  //neigh_bin_linear.resize(no_bins.x * no_bins.y * no_bins.z);

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
        if (domain->boundary_condition.x == 0)
        {
          if (ind_x_min == -1)
            ind_x_min = 0;
          if (ind_x_max == no_bins.x)
            ind_x_max = no_bins.x - 1;
        }

        if (domain->boundary_condition.y == 0)
        {
          if (ind_y_min == -1)
            ind_y_min = 0;
          if (ind_y_max == no_bins.y)
            ind_y_max = no_bins.y - 1;
        }

        if (domain->boundary_condition.z == 0)
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
              caviar::Vector<int> nvec{ind_x_i, ind_y_i, ind_z_i};

#if defined(CAVIAR_WITH_MPI) && !defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#else

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
              
#endif

              neigh_bin[k].emplace_back(nvec);

              //neigh_bin_linear[k].emplace_back(nvec.x + no_bins.x * nvec.y + no_bins.x * no_bins.y * nvec.z);
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

  // for (unsigned int i = 0; i < neigh_bin_linear.size(); ++i)
  // {
  //   remove_duplicates(neigh_bin_linear[i]);
  // }



}

CAVIAR_NAMESPACE_CLOSE
