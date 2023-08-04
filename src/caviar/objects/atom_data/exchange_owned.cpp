
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

bool Atom_data::exchange_owned(long step)
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

  auto &pos = atom_struct_owned.position;

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  const auto me = domain->me;
  int num_local_atoms = atom_struct_owned.id.size();

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
            atom_struct_owned.msd_domain_cross[i].x -= 1;
        }
        while (pos[i].x > x_lupp)
        {
          pos[i].x -= x_width;
          if (msd_process)
            atom_struct_owned.msd_domain_cross[i].x += 1;
        }
      }
      if (bc.y == 1)
      {
        while (pos[i].y < y_llow)
        {
          pos[i].y += y_width;
          if (msd_process)

            atom_struct_owned.msd_domain_cross[i].y -= 1;
        }
        while (pos[i].y > y_lupp)
        {
          pos[i].y -= y_width;
          if (msd_process)

            atom_struct_owned.msd_domain_cross[i].y += 1;
        }
      }
      if (bc.z == 1)
      {
        while (pos[i].z < z_llow)
        {
          pos[i].z += z_width;
          if (msd_process)

            atom_struct_owned.msd_domain_cross[i].z -= 1;
        }
        while (pos[i].z > z_lupp)
        {
          pos[i].z -= z_width;
          if (msd_process)

            atom_struct_owned.msd_domain_cross[i].z += 1;
        }
      }
    }
  }
#elif defined(CAVIAR_WITH_MPI)
  // MPI_Barrier(mpi_comm);

  auto &vel = atom_struct_owned.velocity;
  auto &acc = atom_struct_owned.acceleration;
  auto &id = atom_struct_owned.id;
  auto &type = atom_struct_owned.type;
  auto &msd = atom_struct_owned.msd_domain_cross;

  const auto grid_index_x = domain->grid_index_x;
  const auto grid_index_y = domain->grid_index_y;
  const auto grid_index_z = domain->grid_index_z;

  const auto nprocs_x = domain->nprocs_x;
  const auto nprocs_y = domain->nprocs_y;
  const auto nprocs_z = domain->nprocs_z;

  const auto me = domain->me;

  {
    MPI_Barrier(mpi_comm);
    int local_pos_size = pos.size();
    int global_pos_size = 0;
    MPI_Allreduce(&local_pos_size,
                  &global_pos_size,
                  1, MPI::DOUBLE, MPI_SUM, MPI::COMM_WORLD);
    // // std::cout << step << " : me : " << me << " A local:" << local_pos_size << " global:" << global_pos_size << std::endl;
  }

  const auto &all = domain->all;

  std::vector<int> send_index[3][3][3]; // the index of std::vector<> of the owned
  std::vector<int> send_index_all;      // the index of std::vector<> of the owned, in a 1D vector

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
    recv_mpi_tag[i][j][k] = 26 - m;
    m++;
    FOR_IJK_LOOP_END
  }
  unsigned num_local_atoms = atom_struct_owned.id.size();

  for (unsigned m = 0; m < num_local_atoms; ++m)
  {
    const auto xlc = pos[m].x < x_llow;
    const auto xuc = pos[m].x > x_lupp;
    const auto ylc = pos[m].y < y_llow;
    const auto yuc = pos[m].y > y_lupp;
    const auto zlc = pos[m].z < z_llow;
    const auto zuc = pos[m].z > z_lupp;

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

    // periodic or non-periodic boundary condition
    if (grid_index_x == 0)
      x_val *= bc.x;
    if (grid_index_x == nprocs_x - 1)
      x_val *= bc.x;
    if (grid_index_y == 0)
      y_val *= bc.y;
    if (grid_index_y == nprocs_y - 1)
      y_val *= bc.y;
    if (grid_index_z == 0)
      z_val *= bc.z;
    if (grid_index_z == nprocs_z - 1)
      z_val *= bc.z;

    if (x_val == 0 && y_val == 0 && z_val == 0)
      continue;
    int i = x_val + 1;
    int j = y_val + 1;
    int k = z_val + 1;

    send_index[i][j][k].emplace_back(i);

    if (me != all[i][j][k])
    {
      send_index_all.emplace_back(i);
      send_num[i][j][k]++;
    }
  }

  // // std::cout << step << " : me : " << me << " X 1" << std::endl;

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
  // N x atom_struct_owned.velocity         [5N    : 8N-1 ]
  // N x atom_struct_owned.acceleration     [8N    : 11N-1]
  // (if msd_process==true):step
  // N x atom_struct_owned.msd_domain_cross [11N   : 14N-1 ]
  // BOND ? ANGLE? DIHEDRAL ?

  std::vector<double> send_data[3][3][3], recv_data[3][3][3];

  auto &mpinf = mpi_packet_info_owned;

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

    mpinf.vel = mpinf.total;
    mpinf.total += 3;

    mpinf.acc = mpinf.total;
    mpinf.total += 3;

    if (record_owned_position_old)
    {
      mpinf.pos_o = mpinf.total;
      mpinf.total += 3;
    }
    if (record_owned_velocity_old)
    {
      mpinf.vel_o = mpinf.total;
      mpinf.total += 3;
    }
    if (record_owned_acceleration_old)
    {
      mpinf.acc_o = mpinf.total;
      mpinf.total += 3;
    }
    if (msd_process)
    {
      mpinf.msd = mpinf.total;
      mpinf.total += 3;
    }

    mpinf.mol_ind = mpinf.total;
    mpinf.total += 1;

    mpinf.atomic_bc = mpinf.total;
    mpinf.total += 1;

    // std::cout << "o: mpinf.id:" << mpinf.id << std::endl;
    // std::cout << "o: mpinf.type:" << mpinf.type << std::endl;
    // std::cout << "o: mpinf.pos:" << mpinf.pos << std::endl;
    // std::cout << "o: mpinf.vel:" << mpinf.vel << std::endl;
    // std::cout << "o: mpinf.acc:" << mpinf.acc << std::endl;
    // std::cout << "o: mpinf.pos_o:" << mpinf.pos_o << std::endl;
    // std::cout << "o: mpinf.vel_o:" << mpinf.vel_o << std::endl;
    // std::cout << "o: mpinf.acc_o:" << mpinf.acc_o << std::endl;
    // std::cout << "o: mpinf.msd:" << mpinf.msd << std::endl;
    // std::cout << "o: mpinf.mol_ind:" << mpinf.mol_ind << std::endl;
    // std::cout << "o: mpinf.atomic_bc:" << mpinf.atomic_bc << std::endl;
    // std::cout << "o: mpinf.total:" << mpinf.total << std::endl;
  }

  // // std::cout << step << " : me : " << me << " X 2" << std::endl;

  // ================================================
  // resize send_data to the correct value to increase filling performance
  // ================================================
  // std::cout << step << " , owned me:" << me << " y 1" << std::endl;

  FOR_IJK_LOOP_START
  send_data[i][j][k].resize(mpinf.total * send_num[i][j][k], 0);
  FOR_IJK_LOOP_END

  // // std::cout << step << " : me : " << me << " X 3" << std::endl;

  // std::cout << step << " , owned me:" << me << " y 2" << std::endl;

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
    send_data[i][j][k][mpinf.id * N + c] = id[m];

    send_data[i][j][k][mpinf.type * N + c] = type[m];

    send_data[i][j][k][(mpinf.pos * N) + (3 * c) + 0] = pos[m].x;
    send_data[i][j][k][(mpinf.pos * N) + (3 * c) + 1] = pos[m].y;
    send_data[i][j][k][(mpinf.pos * N) + (3 * c) + 2] = pos[m].z;

    send_data[i][j][k][(mpinf.vel * N) + (3 * c) + 0] = vel[m].x;
    send_data[i][j][k][(mpinf.vel * N) + (3 * c) + 1] = vel[m].y;
    send_data[i][j][k][(mpinf.vel * N) + (3 * c) + 2] = vel[m].z;

    send_data[i][j][k][(mpinf.acc * N) + (3 * c) + 0] = acc[m].x;
    send_data[i][j][k][(mpinf.acc * N) + (3 * c) + 1] = acc[m].y;
    send_data[i][j][k][(mpinf.acc * N) + (3 * c) + 2] = acc[m].z;

    if (msd_process)
    {
      send_data[i][j][k][(mpinf.msd * N) + (3 * c) + 0] = msd[m].x;
      send_data[i][j][k][(mpinf.msd * N) + (3 * c) + 1] = msd[m].y;
      send_data[i][j][k][(mpinf.msd * N) + (3 * c) + 2] = msd[m].z;
    }

    send_data[i][j][k][(mpinf.mol_ind * N) + (3 * c) + 2] = atom_struct_owned.molecule_index[m];
    send_data[i][j][k][(mpinf.atomic_bc * N) + (3 * c) + 2] = atom_struct_owned.atomic_bond_count[m];

    c++;
  }

  FOR_IJK_LOOP_END

  // // std::cout << step << " : me : " << me << " X 4" << std::endl;
  // std::cout << step << " , owned me:" << me << " y 3" << std::endl;

  // ================================================
  // send num of owned
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  MPI_Send(&send_num[i][j][k], 1, MPI_INT, all[i][j][k], send_mpi_tag[i][j][k], mpi_comm); // TAG 0

  FOR_IJK_LOOP_END

  // // std::cout << step << " : me : " << me << " X 5" << std::endl;

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  MPI_Recv(&recv_num[i][j][k], 1, MPI_INT, all[i][j][k], recv_mpi_tag[i][j][k], mpi_comm, MPI_STATUS_IGNORE); // TAG 0

  FOR_IJK_LOOP_END

  // // std::cout << step << " : me : " << me << " X 6" << std::endl;
  // std::cout << step << " , owned me:" << me << " y 4" << std::endl;

  // ================================================
  // resize recv_data to the correct value to increase filling performance
  // ================================================

  FOR_IJK_LOOP_START
  recv_data[i][j][k].resize(mpinf.total * recv_num[i][j][k], 0);
  FOR_IJK_LOOP_END

  // // std::cout << step << " : me : " << me << " X 7" << std::endl;
  // std::cout << step << " , owned me:" << me << " y 5" << std::endl;

  // ================================================
  // SEND OWNED DATA
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  MPI_Send(send_data[i][j][k].data(), mpinf.total * send_num[i][j][k], MPI_DOUBLE, all[i][j][k], send_mpi_tag[i][j][k], mpi_comm); // TAG 1

  FOR_IJK_LOOP_END

  // // std::cout << step << " : me : " << me << " X 8" << std::endl;

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  MPI_Recv(recv_data[i][j][k].data(), mpinf.total * recv_num[i][j][k], MPI_DOUBLE, all[i][j][k], recv_mpi_tag[i][j][k], mpi_comm, MPI_STATUS_IGNORE); // TAG 1

  FOR_IJK_LOOP_END

  // // std::cout << step << " : me : " << me << " X 9" << std::endl;
  // std::cout << step << " , owned me:" << me << " y 6" << std::endl;

  // ================================================
  // FILL the atom_data.owned with depacketing recv_data
  // ================================================

  FOR_IJK_LOOP_START

  if (me == all[i][j][k])
    continue;
  // std::cout << step << " , owned me:" << me << " y 6.1" << std::endl;

  auto N = recv_num[i][j][k];
  // std::cout << step << " , owned me:" << me << " y 6.2" << std::endl;

  if (N == 0)
    continue;
  // std::cout << step << " , owned me:" << me << " y 6.3" << std::endl;

  auto old_size = id.size();
  // std::cout << step << " , owned me:" << me << " y 6.4" << std::endl;

  auto new_size = old_size + N;

  atom_struct_owned_resize(new_size);
  // std::cout << step << " , owned me:" << me << " y 6.5" << std::endl;

  for (int c = 0; c < N; ++c)
  {
    // std::cout << step << " , owned me:" << me << " y 6.6" << std::endl;

    int m = old_size + c;
    // std::cout << step << " , owned me:" << me << " y 6.7" << std::endl;

    id[m] = recv_data[i][j][k][mpinf.id * N + c];
    // std::cout << step << " , owned me:" << me << " y 6.8" << std::endl;

    type[m] = recv_data[i][j][k][mpinf.type * N + c];
    // std::cout << step << " , owned me:" << me << " y 6.9" << std::endl;

    pos[m].x = recv_data[i][j][k][(mpinf.pos * N) + (3 * c) + 0];
    pos[m].y = recv_data[i][j][k][(mpinf.pos * N) + (3 * c) + 1];
    pos[m].z = recv_data[i][j][k][(mpinf.pos * N) + (3 * c) + 2];

    // std::cout << step << " , owned me:" << me << " y 6.11" << std::endl;

    vel[m].x = recv_data[i][j][k][(mpinf.vel * N) + (3 * c) + 0];
    vel[m].y = recv_data[i][j][k][(mpinf.vel * N) + (3 * c) + 1];
    vel[m].z = recv_data[i][j][k][(mpinf.vel * N) + (3 * c) + 2];
    // std::cout << step << " , owned me:" << me << " y 6.12" << std::endl;

    acc[m].x = recv_data[i][j][k][(mpinf.acc * N) + (3 * c) + 0];
    acc[m].y = recv_data[i][j][k][(mpinf.acc * N) + (3 * c) + 1];
    acc[m].z = recv_data[i][j][k][(mpinf.acc * N) + (3 * c) + 2];
    // std::cout << step << " , owned me:" << me << " y 6.13" << std::endl;

    if (msd_process)
    {
      // std::cout << step << " , owned me:" << me << " y 6.14" << std::endl;

      msd[m].x = recv_data[i][j][k][(mpinf.msd * N) + (3 * c) + 0];
      msd[m].y = recv_data[i][j][k][(mpinf.msd * N) + (3 * c) + 1];
      msd[m].z = recv_data[i][j][k][(mpinf.msd * N) + (3 * c) + 2];
    }
    // std::cout << step << " , owned me:" << me << " y 6.15" << std::endl;

    atom_struct_owned.molecule_index[m] = recv_data[i][j][k][(mpinf.mol_ind * N) + (3 * c) + 2];
    atom_struct_owned.atomic_bond_count[m] = recv_data[i][j][k][(mpinf.atomic_bc * N) + (3 * c) + 2];
    // std::cout << step << " , owned me:" << me << " y 7" << std::endl;

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

  // // std::cout << step << " : me : " << me << " X 10" << std::endl;
  // std::cout << step << " , owned me:" << me << " y 8" << std::endl;

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

  // // std::cout << step << " : me : " << me << " X 11 : send_index_all.size()" << send_index_all.size() << std::endl;
  // std::cout << step << " , owned me:" << me << " y 9" << std::endl;

  // ================================================
  // Deleting the particles which are send to another domains
  // ================================================

  if (send_index_all.size() > 0)
  {
    remove_atom(send_index_all);
    make_neighlist = true;
  }
  // // std::cout << step << " : me : " << me << " X 12" << std::endl;
  // std::cout << step << " , owned me:" << me << " y 10" << std::endl;

  {
    MPI_Barrier(mpi_comm);
    long local_pos_size = pos.size();
    long global_pos_size = 0;
    MPI_Allreduce(&local_pos_size,
                  &global_pos_size,
                  1, MPI::LONG, MPI_SUM, MPI::COMM_WORLD);

    // // std::cout << step << " : me : " << me << " B local:" << local_pos_size << " global:" << global_pos_size << std::endl;
  }
  // std::cout << step << " , owned me:" << me << " y 11" << std::endl;

#else

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
  int num_local_atoms = atom_struct_owned.id.size();

  for (unsigned int i = 0; i < num_local_atoms; ++i)
  {
    if (bc.x == 1)
    {
      while (pos[i].x < x_llow)
      {
        pos[i].x += x_width;

        if (msd_process)
          atom_struct_owned.msd_domain_cross[i].x -= 1;
      }
      while (pos[i].x > x_lupp)
      {
        pos[i].x -= x_width;

        if (msd_process)
          atom_struct_owned.msd_domain_cross[i].x += 1;
      }
    }
    if (bc.y == 1)
    {
      while (pos[i].y < y_llow)
      {
        pos[i].y += y_width;

        if (msd_process)
          atom_struct_owned.msd_domain_cross[i].y -= 1;
      }
      while (pos[i].y > y_lupp)
      {
        pos[i].y -= y_width;

        if (msd_process)
          atom_struct_owned.msd_domain_cross[i].y += 1;
      }
    }
    if (bc.z == 1)
    {
      while (pos[i].z < z_llow)
      {
        pos[i].z += z_width;

        if (msd_process)
          atom_struct_owned.msd_domain_cross[i].z -= 1;
      }
      while (pos[i].z > z_lupp)
      {
        pos[i].z -= z_width;

        if (msd_process)
          atom_struct_owned.msd_domain_cross[i].z += 1;
      }
    }
  }
#endif
  return make_neighlist;
}

CAVIAR_NAMESPACE_CLOSE
