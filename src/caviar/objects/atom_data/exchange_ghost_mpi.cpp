
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
#include "caviar/objects/domain.h"

#include <algorithm>
#ifdef CAVIAR_WITH_MPI
#include <mpi.h>
#endif

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

void Atom_data::exchange_ghost_mpi(long ) // timestep
{
#if defined(CAVIAR_WITH_MPI)

  atom_struct_ghost.position.clear();
  atom_struct_ghost.velocity.clear();
  atom_struct_ghost.id.clear();
  atom_struct_ghost.type.clear();

  const auto bc = domain->boundary_condition;

  const auto x_llow = domain->lower_local.x + ghost_cutoff;
  const auto x_lupp = domain->upper_local.x - ghost_cutoff;
  const auto y_llow = domain->lower_local.y + ghost_cutoff;
  const auto y_lupp = domain->upper_local.y - ghost_cutoff;
  const auto z_llow = domain->lower_local.z + ghost_cutoff;
  const auto z_lupp = domain->upper_local.z - ghost_cutoff;

  const auto x_width = domain->upper_local.x - domain->lower_local.x;
  const auto y_width = domain->upper_local.y - domain->lower_local.y;
  const auto z_width = domain->upper_local.z - domain->lower_local.z;

  auto &pos = atom_struct_owned.position;
  auto &vel = atom_struct_owned.velocity;
  auto &id = atom_struct_owned.id;
  auto &type = atom_struct_owned.type;

  auto &g_pos = atom_struct_ghost.position;
  auto &g_vel = atom_struct_ghost.velocity;
  auto &g_id = atom_struct_ghost.id;
  auto &g_type = atom_struct_ghost.type;

 

  const auto grid_index_x = domain->grid_index_x;
  const auto grid_index_y = domain->grid_index_y;
  const auto grid_index_z = domain->grid_index_z;

  const auto nprocs_x = domain->nprocs_x;
  const auto nprocs_y = domain->nprocs_y;
  const auto nprocs_z = domain->nprocs_z;

  //  const auto nprocs = domain->nprocs;
  const auto me = domain->me;
  const auto &all = domain->all;


  std::vector<int> send_id[3][3][3];    // global id of the atom: atom_struct_owned.id
  std::vector<int> recv_id[3][3][3];    // global id of the atom: atom_struct_owned.id
  std::vector<int> send_index[3][3][3]; // the index of std::vector<> of the owned

  int send_num[3][3][3]; // num of owned to be send to the domain all[i][j][k]
  int recv_num[3][3][3]; // num of owned to be recieved from domain all[i][j][k]

  int send_mpi_tag[3][3][3]; // since there might be two messages from the same domain to another but from different angles,
  int recv_mpi_tag[3][3][3]; // , this tag helps to distinguish messages form each other.
  {
    int m = 0;
    FOR_IJK_LOOP_START
    send_num[i][j][k] = 0;
    recv_num[i][j][k] = 0;
    send_mpi_tag[i][j][k] = m;
    recv_mpi_tag[i][j][k] = m;
    m++;
    FOR_IJK_LOOP_END
  }

  unsigned pos_size = pos.size();

  for (unsigned i = 0; i < pos_size; ++i)
  {
    

    //if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
    //  continue;


    const auto xlc = pos[i].x < x_llow;
    const auto xuc = pos[i].x > x_lupp;
    const auto ylc = pos[i].y < y_llow;
    const auto yuc = pos[i].y > y_lupp;
    const auto zlc = pos[i].z < z_llow;
    const auto zuc = pos[i].z > z_lupp;

    int x_val = 0, y_val = 0, z_val = 0;
    if (xlc)
    {
      x_val = -1;
      if (grid_index_x == 0)
        x_val *= bc.x;
    }
    else if (xuc)
    {
      x_val = +1;
      if (grid_index_x == nprocs_x - 1)
        x_val *= bc.x;       
    }

    if (ylc)
    {
      y_val = -1;
      if (grid_index_y == 0)
        y_val *= bc.y;       
    }
    else if (yuc)
    {
      y_val = +1;
      if (grid_index_y == nprocs_y - 1)
        y_val *= bc.y;       
    }

    if (zlc)
    {
      z_val = -1;
      if (grid_index_z == 0)
        z_val *= bc.z;       
    }
    else if (zuc)
    {
      z_val = +1;
      if (grid_index_z == nprocs_z - 1)
        z_val *= bc.z;       
    }

    if (x_val == 0 && y_val == 0 && z_val == 0)
      continue;

    if (x_val != 0)
    {
      send_id[x_val + 1][1][1].emplace_back(id[i]);
      send_index[x_val + 1][1][1].emplace_back(i);
      send_num[x_val + 1][1][1]++;
    }
    if (y_val != 0)
    {
      send_id[1][y_val + 1][1].emplace_back(id[i]);
      send_index[1][y_val + 1][1].emplace_back(i);
      send_num[1][y_val + 1][1]++;
    }
    if (z_val != 0)
    {
      send_id[1][1][z_val + 1].emplace_back(id[i]);
      send_index[1][1][z_val + 1].emplace_back(i);
      send_num[1][1][z_val + 1]++;
    }
    if (x_val != 0 && y_val != 0)
    {
      send_id[x_val + 1][y_val + 1][1].emplace_back(id[i]);
      send_index[x_val + 1][y_val + 1][1].emplace_back(i);
      send_num[x_val + 1][y_val + 1][1]++;
    }
    if (x_val != 0 && z_val != 0)
    {
      send_id[x_val + 1][1][z_val + 1].emplace_back(id[i]);
      send_index[x_val + 1][1][z_val + 1].emplace_back(i);
      send_num[x_val + 1][1][z_val + 1]++;
    }
    if (y_val != 0 && z_val != 0)
    {
      send_id[1][y_val + 1][z_val + 1].emplace_back(id[i]);
      send_index[1][y_val + 1][z_val + 1].emplace_back(i);
      send_num[1][y_val + 1][z_val + 1]++;
    }
    if (x_val != 0 && y_val != 0 && z_val != 0)
    {
      send_id[x_val + 1][y_val + 1][z_val + 1].emplace_back(id[i]);
      send_index[x_val + 1][y_val + 1][z_val + 1].emplace_back(i);
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
  // N x atom_struct_owned.id               [0     : N-1  ]  0,1,2,3
  // N x atom_struct_owned.type             [N     : 2N-1 ]  4,5,6,7
  // N x atom_struct_owned.position         [2N    : 5N-1 ] (8,9,10),(11,12,13),(14,15,16),(17,18,19)
  // if (make_ghost_velocity):
  // N x atom_struct_owned.velocity         [5N    : 8N-1 ]
  // BOND ? ANGLE? DIHEDRAL ?

  std::vector<double> send_data[3][3][3], recv_data[3][3][3];

  auto &mpinf = mpi_packet_info_ghost;

  if (!mpinf.initialized)
  {
    mpinf.initialized = true;
    mpinf.total = 0;

    mpinf.id = mpinf.total;
    mpinf.total += 1;

    mpinf.type = mpinf.total;
    mpinf.total += 1;

    mpinf.pos = mpinf.total;
    mpinf.total += 3;

    if (make_ghost_velocity)
    {
      mpinf.vel = mpinf.total;
      mpinf.total += 3;
    }


  }

  // ================================================
  // resize send_data to the correct value to increase filling performance
  // ================================================

  FOR_IJK_LOOP_START
  if (send_num[i][j][k] > 0)
    send_data[i][j][k].resize(mpinf.total * send_num[i][j][k], 0);
  FOR_IJK_LOOP_END

  // ================================================
  //  FILL the send_data packets
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  int c = 0;
  auto N = send_num[i][j][k];
  if (N == 0)
    continue;
  for (auto m : send_index[i][j][k])
  {
    send_data[i][j][k][mpinf.id * N + c] = id[m];

    send_data[i][j][k][mpinf.type * N + c] = type[m];

    send_data[i][j][k][(mpinf.pos * N) + (3 * c) + 0] = pos[m].x;
    send_data[i][j][k][(mpinf.pos * N) + (3 * c) + 1] = pos[m].y;
    send_data[i][j][k][(mpinf.pos * N) + (3 * c) + 2] = pos[m].z;

    if (make_ghost_velocity)
    {
      send_data[i][j][k][(mpinf.vel * N) + (3 * c) + 0] = vel[m].x;
      send_data[i][j][k][(mpinf.vel * N) + (3 * c) + 1] = vel[m].y;
      send_data[i][j][k][(mpinf.vel * N) + (3 * c) + 2] = vel[m].z;
    }
    c++;
  }

  FOR_IJK_LOOP_END

  // ================================================
  // send num of ghosts
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  MPI_Send(&send_num[i][j][k], 1, MPI_INT, all[i][j][k], send_mpi_tag[i][j][k], MPI_COMM_WORLD); // TAG 0

  FOR_IJK_LOOP_END

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  MPI_Recv(&recv_num[i][j][k], 1, MPI_INT, all[i][j][k], recv_mpi_tag[i][j][k], MPI_COMM_WORLD, MPI_STATUS_IGNORE); // TAG 0

  FOR_IJK_LOOP_END


  // ================================================
  // resize recv_data to the correct value to increase filling performance
  // ================================================

  FOR_IJK_LOOP_START
  if (recv_num[i][j][k] > 0)
  {
    recv_data[i][j][k].resize(mpinf.total * recv_num[i][j][k], 0);
  }


  FOR_IJK_LOOP_END

  // ================================================
  // SEND OWNED DATA
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;
  if (send_num[i][j][k] > 0)
  {

    MPI_Send(send_data[i][j][k].data(), mpinf.total * send_num[i][j][k], MPI_DOUBLE, all[i][j][k], send_mpi_tag[i][j][k], MPI_COMM_WORLD); // TAG 1
  }
  FOR_IJK_LOOP_END

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  if (recv_num[i][j][k] > 0)
  {
    MPI_Recv(recv_data[i][j][k].data(), mpinf.total * recv_num[i][j][k], MPI_DOUBLE, all[i][j][k], recv_mpi_tag[i][j][k], MPI_COMM_WORLD, MPI_STATUS_IGNORE); // TAG 1

  }
  FOR_IJK_LOOP_END

  // ================================================
  // FILL the atom_data.owned with depacketing recv_data
  // ================================================

  FOR_IJK_LOOP_START

  if (me == all[i][j][k])
    continue;

  auto N = recv_num[i][j][k];

  if (N == 0)
    continue;

  auto old_size = g_pos.size(); // g_id.size();
  auto new_size = old_size + N;
  if (N == 0)
    continue;
  atom_struct_ghost_resize(new_size);

  g_id.resize(new_size);
  g_type.resize(new_size);
  g_pos.resize(new_size);
  if (make_ghost_velocity)
    g_vel.resize(new_size);


  for (int c = 0; c < N; ++c)
  {

    int m = old_size + c;


    g_id[m] = recv_data[i][j][k][mpinf.id * N + c];

    g_type[m] = recv_data[i][j][k][mpinf.type * N + c];

    g_pos[m].x = recv_data[i][j][k][(mpinf.pos * N) + (3 * c) + 0];
    g_pos[m].y = recv_data[i][j][k][(mpinf.pos * N) + (3 * c) + 1];
    g_pos[m].z = recv_data[i][j][k][(mpinf.pos * N) + (3 * c) + 2];

    if (make_ghost_velocity)
    {
      g_vel[m].x = recv_data[i][j][k][(mpinf.vel * N) + (3 * c) + 0];
      g_vel[m].y = recv_data[i][j][k][(mpinf.vel * N) + (3 * c) + 1];
      g_vel[m].z = recv_data[i][j][k][(mpinf.vel * N) + (3 * c) + 2];
    }

    // ================================================
    // Applying periodic boundary condition for particles comming from other domains
    // ================================================

    if (i == 0)
      while (g_pos[m].x < x_llow)
        g_pos[m].x += x_width;
    if (j == 0)
      while (g_pos[m].y < y_llow)
        g_pos[m].y += y_width;
    if (k == 0)
      while (g_pos[m].z < z_llow)
        g_pos[m].z += z_width;
    if (i == 2)
      while (g_pos[m].x > x_lupp)
        g_pos[m].x -= x_width;
    if (j == 2)
      while (g_pos[m].y > y_lupp)
        g_pos[m].y -= y_width;
    if (k == 2)
      while (g_pos[m].z > z_lupp)
        g_pos[m].z -= z_width;
  }

  FOR_IJK_LOOP_END

  // ================================================
  // Applying periodic boundary condition for ghosts comming from the current domain itself
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
  {
    auto x_val = i - 1, y_val = j - 1, z_val = k - 1;
    for (auto m : send_index[i][j][k])
    {
      g_vel.emplace_back(vel[m].x, vel[m].y, vel[m].z);
      g_pos.emplace_back(pos[m].x - x_val * x_width, pos[m].y - y_val * y_width, pos[m].z - z_val * z_width);
      g_id.emplace_back(id[m]);
      g_type.emplace_back(type[m]);
    }
  }

  //std::cout << "g_pos size " << g_pos.size() << std::endl;

  FOR_IJK_LOOP_END


// MPI_Barrier(MPI_COMM_WORLD);


#endif
}


CAVIAR_NAMESPACE_CLOSE
