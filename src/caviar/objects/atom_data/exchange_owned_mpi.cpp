
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


//======================================================
//                                                    ||
//                                                    ||
//======================================================

bool Atom_data::exchange_owned_mpi(long) // timestep
{
  #if defined(CAVIAR_WITH_MPI)

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
  // MPI_Barrier(mpi_comm);

  auto &vel = atom_struct_owned.velocity;
  auto &acc = atom_struct_owned.acceleration;
  auto &id = atom_struct_owned.id;
  auto &type = atom_struct_owned.type;
  auto &msd = atom_struct_owned.msd_domain_cross;
  auto &molecule_index = atom_struct_owned.molecule_index;

  const auto grid_index_x = domain->grid_index_x;
  const auto grid_index_y = domain->grid_index_y;
  const auto grid_index_z = domain->grid_index_z;

  const auto nprocs_x = domain->nprocs_x;
  const auto nprocs_y = domain->nprocs_y;
  const auto nprocs_z = domain->nprocs_z;

  const auto me = domain->me;

  const auto &all = domain->all;

  std::vector<int> send_index[3][3][3]; // the index of std::vector<> of the owned
  std::vector<int> send_index_all;      // the index of std::vector<> of the owned, in a 1D vector

  //std::vector<int> send_index_molecule[3][3][3];           // the index of std::vector<> of the owned molecules
  std::vector<int> send_index_molecule_send_data[3][3][3]; // the index of std::vector<> of the owned molecules
  std::vector<int> send_index_molecule_recv_data[3][3][3]; // the index of std::vector<> of the owned molecules

  int send_num[3][3][3];           // num of owned to be send to the domain all[i][j][k]
  //int send_num_molecules[3][3][3]; // num of molecules to be send to the domain all[i][j][k]

  int recv_num[3][3][3]; // num of owned to be recieved from domain all[i][j][k]/ used in he

  int send_mpi_tag[3][3][3]; // since there might be two messages from the same domain to another but from different angles,
  int recv_mpi_tag[3][3][3]; // , this tag helps to distinguish messages form each other.
  {
    int m = 0;
    FOR_IJK_LOOP_START
    send_num[i][j][k] = 0;
    recv_num[i][j][k] = 0;
    //send_num_molecules[i][j][k] = 0;
    send_mpi_tag[i][j][k] = m;
    recv_mpi_tag[i][j][k] = 26 - m;
    m++;
    FOR_IJK_LOOP_END
  }

  // ================================================
  // finding the atoms to be send to other domains. These atoms are not part of any molecules.
  // ================================================

  unsigned pos_size = pos.size();

  for (unsigned m = 0; m < pos_size; ++m)
  {
    if (molecule_index[m] > -1)
      continue; // excluding molecules

    if (atom_struct_owned.mpi_rank[m] != my_mpi_rank)
      continue;

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

    send_index[i][j][k].emplace_back(m);

    if (me != all[i][j][k])
    {
      send_index_all.emplace_back(m);
      send_num[i][j][k]++;
    }
  }

  // ================================================
  // finding the molecules to be send to other domains
  // ================================================

  unsigned num_local_molecules = molecule_struct_owned.size();

  for (unsigned m = 0; m < num_local_molecules; ++m)
  {

    if (molecule_struct_owned[m].ghost) continue;

    Vector<double> cm{0, 0, 0}; // center of molecule
    int molecule_size = molecule_struct_owned[m].atom_list.size();
    for (int n = 0; n < molecule_size; ++n)
    {
      int atom_id = molecule_struct_owned[m].atom_list[n];
      cm += pos[atom_id_to_index[atom_id]];
    }
    cm /= molecule_size;

    const auto xlc = cm.x < x_llow;
    const auto xuc = cm.x > x_lupp;
    const auto ylc = cm.y < y_llow;
    const auto yuc = cm.y > y_lupp;
    const auto zlc = cm.z < z_llow;
    const auto zuc = cm.z > z_lupp;

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

    for (int n = 0; n < molecule_size; ++n)
    {
      int atom_id = molecule_struct_owned[m].atom_list[n];
      int ix = atom_id_to_index[atom_id];

      send_index[i][j][k].emplace_back(ix);

      if (me != all[i][j][k])
      {
        send_index_all.emplace_back(ix);
        send_num[i][j][k]++;

        // ======================================
        // Making the packet of molecule data
        // ======================================
        continue; // XXXXXXXXXXXXX
        int  atom_list_size = 0;//atom_list.size(); XXXXXXXXXXXX
        int  atomic_bond_vector_size = molecule_struct_owned[m].atomic_bond_vector.size();
        int  atomic_angle_vector_size = molecule_struct_owned[m].atomic_angle_vector.size();
        int  atomic_properdihedral_vector_size = molecule_struct_owned[m].atomic_properdihedral_vector.size();

        int m_size = send_index_molecule_send_data[i][j][k].size();
        int p = m_size;
        int new_size = m_size + 4 + atom_list_size + (atomic_bond_vector_size * 4) + (atomic_angle_vector_size * 5) + (atomic_properdihedral_vector_size * 5);

        send_index_molecule_send_data[i][j][k].resize(new_size);

        send_index_molecule_send_data[i][j][k][p++] = atom_list_size;
        send_index_molecule_send_data[i][j][k][p++] = atomic_bond_vector_size;
        send_index_molecule_send_data[i][j][k][p++] = atomic_angle_vector_size;
        send_index_molecule_send_data[i][j][k][p++] = atomic_properdihedral_vector_size;

        for (int o = 0; o < atom_list_size; ++o)
        {
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atom_list[o];
        }

        for (int o = 0; o < atomic_bond_vector_size; ++o)
        {
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_bond_vector[o].id_1;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_bond_vector[o].id_2;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_bond_vector[o].type;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_bond_vector[o].length;
        }

        for (int o = 0; o < atomic_angle_vector_size; ++o)
        {
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_angle_vector[o].id_1;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_angle_vector[o].id_2;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_angle_vector[o].id_3;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_angle_vector[o].type;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_angle_vector[o].value;
        }

        for (int o = 0; o < atomic_properdihedral_vector_size; ++o)
        {
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_properdihedral_vector[o].id_1;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_properdihedral_vector[o].id_2;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_properdihedral_vector[o].id_3;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_properdihedral_vector[o].id_4;
          send_index_molecule_send_data[i][j][k][p++] = molecule_struct_owned[m].atomic_properdihedral_vector[o].type;
        }
      }
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


    int mi = atom_struct_owned.molecule_index[m];
    if (mi  > -1) molecule_struct_owned[mi].ghost = true;

    c++;
  }

  FOR_IJK_LOOP_END

  // ================================================
  // send num of owned
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  MPI_Send(&send_num[i][j][k], 1, MPI_INT, all[i][j][k], send_mpi_tag[i][j][k], mpi_comm); // TAG 0

  FOR_IJK_LOOP_END

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  MPI_Recv(&recv_num[i][j][k], 1, MPI_INT, all[i][j][k], recv_mpi_tag[i][j][k], mpi_comm, MPI_STATUS_IGNORE); // TAG 0

  FOR_IJK_LOOP_END

  // ================================================
  // resize recv_data to the correct value to increase filling performance
  // ================================================

  FOR_IJK_LOOP_START
  if (recv_num[i][j][k] > 0)
    recv_data[i][j][k].resize(mpinf.total * recv_num[i][j][k], 0);
  FOR_IJK_LOOP_END

  // ================================================
  // SEND OWNED DATA
  // ================================================

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;
  if (send_num[i][j][k] > 0)
    MPI_Send(send_data[i][j][k].data(), mpinf.total * send_num[i][j][k], MPI_DOUBLE, all[i][j][k], send_mpi_tag[i][j][k], mpi_comm); // TAG 1

  FOR_IJK_LOOP_END

  FOR_IJK_LOOP_START
  if (me == all[i][j][k])
    continue;

  if (recv_num[i][j][k] > 0)
    MPI_Recv(recv_data[i][j][k].data(), mpinf.total * recv_num[i][j][k], MPI_DOUBLE, all[i][j][k], recv_mpi_tag[i][j][k], mpi_comm, MPI_STATUS_IGNORE); // TAG 1

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


  auto old_size = id.size();
  auto new_size = old_size + N;
  atom_struct_owned_resize(new_size);

  for (int c = 0; c < N; ++c)
  {


    int m = old_size + c;
    id[m] = recv_data[i][j][k][mpinf.id * N + c];
    type[m] = recv_data[i][j][k][mpinf.type * N + c];


    pos[m].x = recv_data[i][j][k][(mpinf.pos * N) + (3 * c) + 0];
    pos[m].y = recv_data[i][j][k][(mpinf.pos * N) + (3 * c) + 1];
    pos[m].z = recv_data[i][j][k][(mpinf.pos * N) + (3 * c) + 2];

    vel[m].x = recv_data[i][j][k][(mpinf.vel * N) + (3 * c) + 0];
    vel[m].y = recv_data[i][j][k][(mpinf.vel * N) + (3 * c) + 1];
    vel[m].z = recv_data[i][j][k][(mpinf.vel * N) + (3 * c) + 2];

    acc[m].x = recv_data[i][j][k][(mpinf.acc * N) + (3 * c) + 0];
    acc[m].y = recv_data[i][j][k][(mpinf.acc * N) + (3 * c) + 1];
    acc[m].z = recv_data[i][j][k][(mpinf.acc * N) + (3 * c) + 2];

    if (msd_process)
    {

      msd[m].x = recv_data[i][j][k][(mpinf.msd * N) + (3 * c) + 0];
      msd[m].y = recv_data[i][j][k][(mpinf.msd * N) + (3 * c) + 1];
      msd[m].z = recv_data[i][j][k][(mpinf.msd * N) + (3 * c) + 2];
    }

    atom_struct_owned.molecule_index[m] = recv_data[i][j][k][(mpinf.mol_ind * N) + (3 * c) + 2];
    atom_struct_owned.atomic_bond_count[m] = recv_data[i][j][k][(mpinf.atomic_bc * N) + (3 * c) + 2];

    int mi = atom_struct_owned.molecule_index[m];
    if (mi  > -1) molecule_struct_owned[mi].ghost = false;

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

  //==========================================
  // Checking num_of_particles code snippet
  //==========================================
  // {
  // MPI_Barrier(mpi_comm);
  // long local_pos_size = pos.size();
  // long global_pos_size = 0;
  // MPI_Allreduce(&local_pos_size,
  //               &global_pos_size,
  //               1, MPI::LONG, MPI_SUM, MPI::COMM_WORLD);
  // }

  return make_neighlist;
#else
  return false;
#endif
}



CAVIAR_NAMESPACE_CLOSE
