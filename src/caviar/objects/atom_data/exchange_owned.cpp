
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

#include "caviar/objects/atom_data.h"
#include "caviar/interpreter/communicator.h"
#include "caviar/interpreter/error.h"
#include "caviar/objects/domain.h"

#include <algorithm>

CAVIAR_NAMESPACE_OPEN

#define FOR_IJK_LOOP_START         \
  for (auto i = 0; i < 3; ++i)     \
  {                                \
    for (auto j = 0; j < 3; ++j)   \
    {                              \
      for (auto k = 0; k < 3; ++k) \
      {

#define FOR_IJK_LOOP_END \
  }                      \
  }                      \
  }

bool Atom_data::exchange_owned()
{
  if (domain == nullptr)
    error->all("Atom_data::exchange_owned: domain = nullptr");
  bool make_neighlist = false;

  const auto bc = domain->boundary_condition;

  // cutoff extra may make problem for induced_charge mesh.
  // in the situations that the mesh size is exactly equal to the domain
  // measurements, cutoff_extra may makes the particles go outside of the mesh.
  // that throws an exception when calculating the force or accelerations on
  // the particle. in order to avoid that, one can fix it by setting cutoff_extra
  // to zero, or by fixing particle position in the forcefield calculations.
  const auto x_llow = domain->lower_local.x - cutoff_extra;
  const auto x_lupp = domain->upper_local.x + cutoff_extra;
  const auto y_llow = domain->lower_local.y - cutoff_extra;
  const auto y_lupp = domain->upper_local.y + cutoff_extra;
  const auto z_llow = domain->lower_local.z - cutoff_extra;
  const auto z_lupp = domain->upper_local.z + cutoff_extra;

  const auto x_width = domain->upper_local.x - domain->lower_local.x;
  const auto y_width = domain->upper_local.y - domain->lower_local.y;
  const auto z_width = domain->upper_local.z - domain->lower_local.z;

  auto &pos = owned.position;

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  const auto me = domain->me;
  if (me == 0)
  {
    for (unsigned int i = 0; i < num_local_atoms; ++i)
    {
      if (bc.x == 1)
      { // TODO While should be changed to if in the future
        while (pos[i].x < x_llow)
        {
          pos[i].x += x_width;
          if (msd_process)
            owned.msd_domain_cross[i].x -= 1;
        }
        while (pos[i].x > x_lupp)
        {
          pos[i].x -= x_width;
          if (msd_process)
            owned.msd_domain_cross[i].x += 1;
        }
      }
      if (bc.y == 1)
      {
        while (pos[i].y < y_llow)
        {
          pos[i].y += y_width;
          if (msd_process)

            owned.msd_domain_cross[i].y -= 1;
        }
        while (pos[i].y > y_lupp)
        {
          pos[i].y -= y_width;
          if (msd_process)

            owned.msd_domain_cross[i].y += 1;
        }
      }
      if (bc.z == 1)
      {
        while (pos[i].z < z_llow)
        {
          pos[i].z += z_width;
          if (msd_process)

            owned.msd_domain_cross[i].z -= 1;
        }
        while (pos[i].z > z_lupp)
        {
          pos[i].z -= z_width;
          if (msd_process)

            owned.msd_domain_cross[i].z += 1;
        }
      }
    }
  }
#elif defined(CAVIAR_WITH_MPI)
  //MPI_Barrier(mpi_comm);


  auto &vel = owned.velocity;
  auto &acc = owned.acceleration;
  auto &id = owned.id;
  auto &type = owned.type;
  auto &msd = owned.msd_domain_cross;

  const auto grid_index_x = domain->grid_index_x;
  const auto grid_index_y = domain->grid_index_y;
  const auto grid_index_z = domain->grid_index_z;

  const auto nprocs_x = domain->nprocs_x;
  const auto nprocs_y = domain->nprocs_y;
  const auto nprocs_z = domain->nprocs_z;

  const auto me = domain->me;


  const auto &all = domain->all;

  std::vector<int> send_id[3][3][3];    // global id of the atom: owned.id
  std::vector<int> recv_id[3][3][3];    // global id of the atom: owned.id
  std::vector<int> send_index[3][3][3]; // the index of std::vector<> of the owned
  std::vector<int> send_index_all;      // the index of std::vector<> of the owned, in a 1D vector

  int send_num[3][3][3]; // num of owned to be send to the domain all[i][j][k]
  int recv_num[3][3][3]; // num of owned to be recieved from domain all[i][j][k]


  FOR_IJK_LOOP_START
  send_num[i][j][k] = 0;
  recv_num[i][j][k] = 0;
  FOR_IJK_LOOP_END


  for (unsigned i = 0; i < num_local_atoms; ++i)
  {
    const auto xlc = pos[i].x < x_llow;
    const auto xuc = pos[i].x > x_lupp;
    const auto ylc = pos[i].y < y_llow;
    const auto yuc = pos[i].y > y_lupp;
    const auto zlc = pos[i].z < z_llow;
    const auto zuc = pos[i].z > z_lupp;

    int x_val = 0, y_val = 0, z_val = 0;

    if (xlc)
      x_val = -1;
    if (xuc)
      x_val = +1;
    if (ylc)
      y_val = -1;
    if (yuc)
      y_val = +1;
    if (zlc)
      z_val = -1;
    if (zuc)
      z_val = +1;

    if (grid_index_x == 0)
      x_val *= bc.x; // periodic or non-periodic boundary condition
    if (grid_index_x == nprocs_x - 1)
      x_val *= bc.x; // //
    if (grid_index_y == 0)
      y_val *= bc.y; // //
    if (grid_index_y == nprocs_y - 1)
      y_val *= bc.y; // //
    if (grid_index_z == 0)
      z_val *= bc.z; // //
    if (grid_index_z == nprocs_z - 1)
      z_val *= bc.z; // //

    if (x_val == 0 && y_val == 0 && z_val == 0)
    {
      continue;
    }
    else
    {
      send_id[x_val + 1][y_val + 1][z_val + 1].emplace_back(id[i]);
      send_index[x_val + 1][y_val + 1][z_val + 1].emplace_back(i);
      send_index_all.emplace_back(i);
      send_num[x_val + 1][y_val + 1][z_val + 1]++;
    }
  }

  // ================================================
  // making the send_data
  // ================================================

  // o_send_data:  contains all the data of the atoms that is transfered to another MPI process
  // All of the data are going to be casted as 'double' type
  // N: the number of atoms that are going to be send to another process
  //
  // The packet is as follows with N=4 particles example
  // N x owned.id               [0     : N-1  ]  0,1,2,3
  // N x owned.type             [N     : 2N-1 ]  4,5,6,7
  // N x owned.position         [2N    : 5N-1 ] (8,9,10),(11,12,13),(14,15,16),(17,18,19)
  // N x owned.velocity         [5N    : 8N-1 ]
  // N x owned.acceleration     [8N    : 11N-1]
  // (if msd_process==true):
  // N x owned.msd_domain_cross [11N   : 14N-1 ]
  // BOND ? ANGLE? DIHEDRAL ?

  std::vector<double> send_data[3][3][3], recv_data[3][3][3];

  // ================================================
  // resize send_data to the correct value to increase filling performance
  // ================================================

  FOR_IJK_LOOP_START
  int N_tot;
  if (msd_process)
    N_tot = 14 * send_num[i][j][k];
  else
    N_tot = 11 * send_num[i][j][k];

  send_data[i][j][k].resize(N_tot, 0);

  FOR_IJK_LOOP_END

  // ================================================
  //  FILL the send_data packets
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  int c = 0;
  auto N = send_num[i][j][k];
  for (auto m : send_index[i][j][k])
  {
    send_data[i][j][k][c] = id[m];

    send_data[i][j][k][N + c] = type[m];

    send_data[i][j][k][(2 * N) + (3 * c) + 0] = pos[m].x;
    send_data[i][j][k][(2 * N) + (3 * c) + 1] = pos[m].y;
    send_data[i][j][k][(2 * N) + (3 * c) + 2] = pos[m].z;

    send_data[i][j][k][(5 * N) + (3 * c) + 0] = vel[m].x;
    send_data[i][j][k][(5 * N) + (3 * c) + 1] = vel[m].y;
    send_data[i][j][k][(5 * N) + (3 * c) + 2] = vel[m].z;

    send_data[i][j][k][(8 * N) + (3 * c) + 0] = acc[m].x;
    send_data[i][j][k][(8 * N) + (3 * c) + 1] = acc[m].y;
    send_data[i][j][k][(8 * N) + (3 * c) + 2] = acc[m].z;

    if (msd_process)
    {
      send_data[i][j][k][(11 * N) + (3 * c) + 0] = msd[m].x;
      send_data[i][j][k][(11 * N) + (3 * c) + 1] = msd[m].y;
      send_data[i][j][k][(11 * N) + (3 * c) + 2] = msd[m].z;
    }

    c++;
  }

  FOR_IJK_LOOP_END

  // ================================================
  // send num of owned
  // ================================================

  //MPI_Barrier(mpi_comm);

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  MPI_Send(&send_num[i][j][k], 1, MPI_INT, all[i][j][k], 0, mpi_comm); // TAG 0

  FOR_IJK_LOOP_END

  //MPI_Barrier(mpi_comm);

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  MPI_Recv(&recv_num[i][j][k], 1, MPI_INT, all[i][j][k], 0, mpi_comm, MPI_STATUS_IGNORE); // TAG 0

  FOR_IJK_LOOP_END


  // ================================================

  // MPI_Barrier(mpi_comm);

  // ================================================
  // resize recv_data to the correct value to increase filling performance
  // ================================================

  FOR_IJK_LOOP_START
  int N_tot;
  if (msd_process)
    N_tot = 14 * recv_num[i][j][k];
  else
    N_tot = 11 * recv_num[i][j][k];

  if (N_tot == 0)
    continue; //recv_data[i][j][k].resize(1, 0);
  else
    recv_data[i][j][k].resize(N_tot, 0);

  FOR_IJK_LOOP_END


  // ================================================
  // SEND OWNED DATA
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  int N_tot;
  if (msd_process)
    N_tot = 14 * send_num[i][j][k];
  else
    N_tot = 11 * send_num[i][j][k];

  MPI_Send(send_data[i][j][k].data(), N_tot, MPI_DOUBLE, all[i][j][k], 1, mpi_comm); // TAG 1

  FOR_IJK_LOOP_END

  MPI_Barrier(mpi_comm);

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  int N_tot;
  if (msd_process)
    N_tot = 14 * recv_num[i][j][k];
  else
    N_tot = 11 * recv_num[i][j][k];

  MPI_Recv(recv_data[i][j][k].data(), N_tot, MPI_DOUBLE, all[i][j][k], 1, mpi_comm, MPI_STATUS_IGNORE); // TAG 1

  FOR_IJK_LOOP_END

  // MPI_Barrier(mpi_comm);

  // ================================================
  // FILL the atom_data.owned with depacketing recv_data
  // ================================================


  FOR_IJK_LOOP_START

  if (me == all[i][j][k])
    continue;

  auto N = recv_num[i][j][k];

  if (N == 0)
    continue;

  auto old_size = id.size();
  auto new_size = old_size + N;
  id.resize(new_size);
  type.resize(new_size);
  pos.resize(new_size);
  vel.resize(new_size);
  acc.resize(new_size);
  if (msd_process)
    msd.resize(new_size);

  num_local_atoms += N;

  for (int c = 0; c < N; ++c)
  {

    int m = old_size + c;

    id[m] = recv_data[i][j][k][c];

    type[m] = recv_data[i][j][k][N + c];

    pos[m].x = recv_data[i][j][k][(2 * N) + (3 * c) + 0];
    pos[m].y = recv_data[i][j][k][(2 * N) + (3 * c) + 1];
    pos[m].z = recv_data[i][j][k][(2 * N) + (3 * c) + 2];

    vel[m].x = recv_data[i][j][k][(5 * N) + (3 * c) + 0];
    vel[m].y = recv_data[i][j][k][(5 * N) + (3 * c) + 1];
    vel[m].z = recv_data[i][j][k][(5 * N) + (3 * c) + 2];

    acc[m].x = recv_data[i][j][k][(8 * N) + (3 * c) + 0];
    acc[m].y = recv_data[i][j][k][(8 * N) + (3 * c) + 1];
    acc[m].z = recv_data[i][j][k][(8 * N) + (3 * c) + 2];

    if (msd_process)
    {
      msd[m].x = recv_data[i][j][k][(11 * N) + (3 * c) + 0];
      msd[m].y = recv_data[i][j][k][(11 * N) + (3 * c) + 1];
      msd[m].z = recv_data[i][j][k][(11 * N) + (3 * c) + 2];
    }


    // ================================================
    // Applying periodic boundary condition for particles comming from other domains
    // ================================================

    if (i == 0)
      while (pos[m].x < x_llow)
        pos[m].x += x_width;
    if (j == 0)
      while (pos[m].y < y_llow)
        pos[m].y += y_width;
    if (k == 0)
      while (pos[m].z < z_llow)
        pos[m].z += z_width;
    if (i == 2)
      while (pos[m].x > x_lupp)
        pos[m].x -= x_width;
    if (j == 2)
      while (pos[m].y > y_lupp)
        pos[m].y -= y_width;
    if (k == 2)
      while (pos[m].z > z_lupp)
        pos[m].z -= z_width;
  }

  FOR_IJK_LOOP_END


  // ================================================
  // Applying periodic boundary condition for particles comming from the current domain itself
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
  {
    auto ii = i - 1, jj = j - 1, kk = k - 1;
    for (auto m : send_index[i][j][k])
    {
      pos[m].x -= ii * x_width;
      pos[m].y -= jj * y_width;
      pos[m].z -= kk * z_width;
    }
  }
  FOR_IJK_LOOP_END


  // ================================================
  // Deleting the particles which are send to another domains
  // ================================================

  if (send_index_all.size() > 0)
  {
    remove_atom(send_index_all);
    make_neighlist = true;
  }
  // MPI_Barrier(mpi_comm);

#else

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < num_local_atoms; ++i)
  {
    if (bc.x == 1)
    {
      while (pos[i].x < x_llow)
      {
        pos[i].x += x_width;
        if (msd_process)

          owned.msd_domain_cross[i].x -= 1;
      }
      while (pos[i].x > x_lupp)
      {
        pos[i].x -= x_width;
        if (msd_process)

          owned.msd_domain_cross[i].x += 1;
      }
    }
    if (bc.y == 1)
    {
      while (pos[i].y < y_llow)
      {
        pos[i].y += y_width;
        if (msd_process)

          owned.msd_domain_cross[i].y -= 1;
      }
      while (pos[i].y > y_lupp)
      {
        pos[i].y -= y_width;
        if (msd_process)

          owned.msd_domain_cross[i].y += 1;
      }
    }
    if (bc.z == 1)
    {
      while (pos[i].z < z_llow)
      {
        pos[i].z += z_width;
        if (msd_process)

          owned.msd_domain_cross[i].z -= 1;
      }
      while (pos[i].z > z_lupp)
      {
        pos[i].z -= z_width;
        if (msd_process)

          owned.msd_domain_cross[i].z += 1;
      }
    }
  }
#endif
  return make_neighlist;
}

CAVIAR_NAMESPACE_CLOSE
