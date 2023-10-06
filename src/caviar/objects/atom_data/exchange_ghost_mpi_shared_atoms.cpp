
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

void Atom_data::exchange_ghost_mpi_shared_atoms(long) // timestep
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

  const auto me = domain->me;
  const auto &all = domain->all;


  auto &send_index = mpi_tools.send_index;  // the index of std::vector<> of the owned
  auto &send_index_all = mpi_tools.send_index_all;      // the index of std::vector<> of the owned, in a 1D vector

  auto &send_data = mpi_tools.send_data; 
  auto &recv_data = mpi_tools.recv_data;
  
  auto &send_num = mpi_tools.send_num; // num of owned to be send to the domain all[i][j][k]
  auto &recv_num = mpi_tools.recv_num; // num of owned to be recieved from domain all[i][j][k]/ used in he

  auto &send_mpi_tag = mpi_tools.send_mpi_tag; // since there might be two messages from the same domain to another but from different angles,
  auto &recv_mpi_tag = mpi_tools.recv_mpi_tag; // , this tag helps to distinguish messages form each other.
  
  MPI_Request mpi_requests[3][3][3];
  
  if (mpi_tools.initialize)
  {
    mpi_tools.initialize = false;
    int m = 0;
    FOR_IJK_LOOP_START
    send_mpi_tag[i][j][k] = m;
    recv_mpi_tag[i][j][k] = 26 - m;
    m++;
    FOR_IJK_LOOP_END
  }

  FOR_IJK_LOOP_START
  send_num[i][j][k] = 0;
  recv_num[i][j][k] = 0;
  send_data[i][j][k].clear();
  recv_data[i][j][k].clear();
  send_index[i][j][k].clear();
  FOR_IJK_LOOP_END

  send_index_all.clear();

  // ================================================
  // finding the atoms to be send to other domains. These atoms are not part of any molecules.
  // ================================================

  unsigned pos_size = pos.size();

  for (unsigned i = 0; i < pos_size; ++i)
  {

    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;

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
      send_index[x_val + 1][1][1].emplace_back(i);
      send_num[x_val + 1][1][1]++;
    }
    if (y_val != 0)
    {
      send_index[1][y_val + 1][1].emplace_back(i);
      send_num[1][y_val + 1][1]++;
    }
    if (z_val != 0)
    {
      send_index[1][1][z_val + 1].emplace_back(i);
      send_num[1][1][z_val + 1]++;
    }
    if (x_val != 0 && y_val != 0)
    {
      send_index[x_val + 1][y_val + 1][1].emplace_back(i);
      send_num[x_val + 1][y_val + 1][1]++;
    }
    if (x_val != 0 && z_val != 0)
    {
      send_index[x_val + 1][1][z_val + 1].emplace_back(i);
      send_num[x_val + 1][1][z_val + 1]++;
    }
    if (y_val != 0 && z_val != 0)
    {
      send_index[1][y_val + 1][z_val + 1].emplace_back(i);
      send_num[1][y_val + 1][z_val + 1]++;
    }
    if (x_val != 0 && y_val != 0 && z_val != 0)
    {
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
  //  FILL the send_data packets
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  int c = 0;
  auto N = send_num[i][j][k];
  if (N == 0)
    continue;

  send_data[i][j][k].resize(mpinf.total * send_num[i][j][k], 0);

  for (auto m : send_index[i][j][k])
  {
    Vector<double> pos_tmp = pos[m];
    //================================================================
    // Fixing the position of the atoms in periodic boundary condition
    //================================================================
    if (bc.x == 1)
    {
      if ((grid_index_x == 0) && (i == 0))
      {
        pos_tmp.x += domain->size_global.x;
      }
    
      if ((grid_index_x == nprocs_x - 1) && (i == 2))
      {
        pos_tmp.x -= domain->size_global.x;
      }
    }

    if (bc.y == 1)
    {
      if ((grid_index_y == 0) && (j == 0))
      {
        pos_tmp.y += domain->size_global.y;
      }

      if ((grid_index_y == nprocs_y - 1) && (j == 2))
      {
        pos_tmp.y -= domain->size_global.y;
      }
    }

    if (bc.z == 1)
    {
      if ((grid_index_z == 0) && (k == 0) )
      {
        pos_tmp.z += domain->size_global.z;
      }

      if ((grid_index_z == nprocs_z - 1) && (k == 2))
      {
        pos_tmp.z -= domain->size_global.z;
      }
    }

    //================================================================
    // Prepairing Send_data
    //================================================================
    send_data[i][j][k][mpinf.id * N + c] = id[m];

    send_data[i][j][k][mpinf.type * N + c] = type[m];

    send_data[i][j][k][(mpinf.pos * N) + (3 * c) + 0] = pos_tmp.x;
    send_data[i][j][k][(mpinf.pos * N) + (3 * c) + 1] = pos_tmp.y;
    send_data[i][j][k][(mpinf.pos * N) + (3 * c) + 2] = pos_tmp.z;

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

  //MPI_Send(&send_num[i][j][k], 1, MPI_INT, all[i][j][k], send_mpi_tag[i][j][k], MPI_COMM_WORLD); // TAG 0
  MPI_Isend(&send_num[i][j][k], 1, MPI_INT, all[i][j][k], send_mpi_tag[i][j][k], MPI_COMM_WORLD, &mpi_requests[i][j][k]); // TAG 0

  FOR_IJK_LOOP_END

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;


  //MPI_Recv(&recv_num[i][j][k], 1, MPI_INT, all[i][j][k], recv_mpi_tag[i][j][k], MPI_COMM_WORLD, MPI_STATUS_IGNORE); // TAG 0
  MPI_Irecv(&recv_num[i][j][k], 1, MPI_INT, all[i][j][k], recv_mpi_tag[i][j][k], MPI_COMM_WORLD, &mpi_requests[i][j][k]); // TAG 0

  FOR_IJK_LOOP_END

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
  continue;
  MPI_Wait(&mpi_requests[i][j][k], MPI_STATUS_IGNORE);
  FOR_IJK_LOOP_END

  // ================================================
  // SEND OWNED DATA
  // ================================================

  FOR_IJK_LOOP_START

  if (me == all[i][j][k])
    continue;

  if (send_num[i][j][k] == 0)
    continue;  

  //MPI_Send(send_data[i][j][k].data(), mpinf.total * send_num[i][j][k], MPI_DOUBLE, all[i][j][k], send_mpi_tag[i][j][k], MPI_COMM_WORLD); // TAG 1
  MPI_Isend(send_data[i][j][k].data(), mpinf.total * send_num[i][j][k], MPI_DOUBLE, all[i][j][k], send_mpi_tag[i][j][k], MPI_COMM_WORLD, &mpi_requests[i][j][k]); // TAG 1
  
  FOR_IJK_LOOP_END

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  if (recv_num[i][j][k] == 0)
    continue;
  
  recv_data[i][j][k].resize(mpinf.total * recv_num[i][j][k], 0);

  //MPI_Recv(recv_data[i][j][k].data(), mpinf.total * recv_num[i][j][k], MPI_DOUBLE, all[i][j][k], recv_mpi_tag[i][j][k], MPI_COMM_WORLD, MPI_STATUS_IGNORE); // TAG 1
  MPI_Irecv(recv_data[i][j][k].data(), mpinf.total * recv_num[i][j][k], MPI_DOUBLE, all[i][j][k], recv_mpi_tag[i][j][k], MPI_COMM_WORLD, &mpi_requests[i][j][k]); // TAG 1
  
  FOR_IJK_LOOP_END

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
  continue;

  if (recv_num[i][j][k] == 0)
  continue;
  MPI_Wait(&mpi_requests[i][j][k], MPI_STATUS_IGNORE);
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
      if (make_ghost_velocity)
        g_vel.emplace_back(vel[m].x, vel[m].y, vel[m].z);
      g_pos.emplace_back(pos[m].x - x_val * domain->size_local.x, pos[m].y - y_val * domain->size_local.y, pos[m].z - z_val * domain->size_local.z);
      g_id.emplace_back(id[m]);
      g_type.emplace_back(type[m]);
    }
  }

  // std::cout << "g_pos size " << g_pos.size() << std::endl;

  FOR_IJK_LOOP_END

  // MPI_Barrier(MPI_COMM_WORLD);


#endif
}

CAVIAR_NAMESPACE_CLOSE
