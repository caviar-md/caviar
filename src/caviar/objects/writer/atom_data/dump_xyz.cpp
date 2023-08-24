
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

#include "caviar/objects/writer/atom_data.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/interpreter/communicator.h"
#include "caviar/utility/time_utility.h"
#include "caviar/objects/unique/time_function_3d.h"
#ifdef CAVIAR_WITH_MPI
#include <mpi.h>
#endif
CAVIAR_NAMESPACE_OPEN

namespace writer
{

  //================================================
  //                                              ||
  //================================================
  void Atom_data::dump_xyz(int64_t i)
  {
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

    if (my_mpi_rank == 0)
    {
      dump_xyz_serial(i);
    }
#elif defined(CAVIAR_WITH_MPI)

    if (comm->nprocs > 1)
    {
      if (mpi_single_file)
      {
         dump_xyz_mpi_shared_atoms(i);
        //dump_xyz_mpi(i);
      }
      if (mpi_separate_files)
        dump_xyz_serial(i, true);
    }
    else
    {
      dump_xyz_serial(i);
    }

#else

    dump_xyz_serial(i);

#endif

    double wallTimeXyzDump2 = get_wall_time();

    double dtstart = wallTimeXyzDump2 - wallTimeXyzDump1;
    std::string s = "writer::atom_data:: dump_xyz at step " + std::to_string(i) +
                    +" . Elapsed time since previous xyz dump: " + std::to_string(dtstart);
    output->info(s, 2);
    wallTimeXyzDump1 = wallTimeXyzDump2;
  }

  //================================================
  //                                              ||
  //================================================
  void Atom_data::dump_xyz_ghost(int64_t i)
  {
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

    if (my_mpi_rank == 0)
    {
      dump_xyz_ghost_serial(i);
    }
#elif defined(CAVIAR_WITH_MPI)

    if (comm->nprocs > 1)
    {
      if (mpi_single_file)
        dump_xyz_ghost_mpi(i);
      if (mpi_separate_files)
        dump_xyz_ghost_serial(i, my_mpi_rank);
    }
    else
    {
      dump_xyz_ghost_serial(i);
    }

#else

    dump_xyz_ghost_serial(i);

#endif
  }
  //================================================
  //                                              ||
  //================================================
  void Atom_data::dump_xyz_mpi(int64_t)
  {
#if defined(CAVIAR_WITH_MPI)

    const auto &pos = atom_data->atom_struct_owned.position;
    const auto &vel = atom_data->atom_struct_owned.velocity;
    const auto &acc = atom_data->atom_struct_owned.acceleration;
    const auto &id = atom_data->atom_struct_owned.id;
    const auto &type = atom_data->atom_struct_owned.type;

    unsigned int send_num = pos.size(); // num_local_atoms;
    const unsigned nprocs = comm->nprocs;
    unsigned int num_total_atoms = pos.size();

    std::vector<unsigned int> recv_num(nprocs, 0);
    //==================================// num_local_atoms send
    if (my_mpi_rank != 0)
    {

      MPI_Send(&send_num, 1, MPI::UNSIGNED, 0, 0, MPI::COMM_WORLD);
    }
    if (my_mpi_rank == 0)
    {
      for (unsigned i = 1; i < nprocs; ++i)
      {
        MPI_Recv(&recv_num[i], 1, MPI::UNSIGNED, i, 0, MPI::COMM_WORLD, MPI_STATUS_IGNORE);

        num_total_atoms += recv_num[i];
      }
    }
    auto &mpinf = mpi_packet_info_owned;

    if (!mpinf.initialized)
    {
      mpinf.initialized = true;
      mpinf.total = 0;
      if (output_id)
      {
        mpinf.id = mpinf.total;
        mpinf.total += 1;
      }
      mpinf.type = mpinf.total;
      mpinf.total += 1;

      mpinf.pos = mpinf.total;
      mpinf.total += 3;

      if (output_velocity)
      {
        mpinf.vel = mpinf.total;
        mpinf.total += 3;
      }
      if (output_acceleration)
      {
        mpinf.acc = mpinf.total;
        mpinf.total += 3;
      }
    }

    std::vector<double> send_data;              // rank != 0
    std::vector<std::vector<double>> recv_data; // rank 0

    if (my_mpi_rank != 0)
    {

      send_data.resize(mpinf.total * send_num, 0);
    }

    if (my_mpi_rank == 0)
    {
      recv_data.resize(nprocs);
      for (unsigned i = 1; i < nprocs; ++i)
      {
        if (recv_num[i] > 0)
          recv_data[i].resize(mpinf.total * recv_num[i], 0);
      }
    }

    if (my_mpi_rank != 0)
    {
      auto N = send_num;

      for (unsigned int i = 0; i < send_num; ++i)
      {
        if (output_id)
          send_data[mpinf.id * N + i] = id[i];

        send_data[mpinf.type * N + i] = type[i];

        send_data[(mpinf.pos * N) + (3 * i) + 0] = pos[i].x;
        send_data[(mpinf.pos * N) + (3 * i) + 1] = pos[i].y;
        send_data[(mpinf.pos * N) + (3 * i) + 2] = pos[i].z;

        if (output_velocity)
        {
          send_data[(mpinf.vel * N) + (3 * i) + 0] = vel[i].x;
          send_data[(mpinf.vel * N) + (3 * i) + 1] = vel[i].y;
          send_data[(mpinf.vel * N) + (3 * i) + 2] = vel[i].z;
        }

        if (output_acceleration)
        {
          send_data[(mpinf.acc * N) + (3 * i) + 0] = acc[i].x;
          send_data[(mpinf.acc * N) + (3 * i) + 1] = acc[i].y;
          send_data[(mpinf.acc * N) + (3 * i) + 2] = acc[i].z;
        }
      }
    }
    //-----------------------------------------------//
    if (my_mpi_rank != 0)
      if (send_num > 0)

        MPI_Send(send_data.data(), mpinf.total * send_num, MPI_DOUBLE, 0, 0, MPI::COMM_WORLD); // TAG 1

    if (my_mpi_rank == 0)
      for (unsigned i = 1; i < nprocs; ++i)
        if (recv_num[i] > 0)

          MPI_Recv(recv_data[i].data(), mpinf.total * recv_num[i], MPI_DOUBLE, i, 0, MPI::COMM_WORLD, MPI_STATUS_IGNORE); // TAG 1

    //-----------------------------------------------//
    if (my_mpi_rank != 0)
      return;

    Vector<double> p_o{0, 0, 0};

    if (position_offset != nullptr)
      p_o = position_offset->current_value;

    ofs_xyz << num_total_atoms << "\nAtom\n";

    for (unsigned int i = 0; i < send_num; ++i)
    {
      ofs_xyz << type[i];
      if (output_id)
      {

        ofs_xyz << " " << id[i];
      }

      ofs_xyz << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;

      if (output_velocity)
        ofs_xyz << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
      if (output_acceleration)
        ofs_xyz << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
      ofs_xyz << "\n";
    }

    for (unsigned i = 1; i < nprocs; ++i)
    {

      int N = recv_num[i];

      for (unsigned int j = 0; j < recv_num[i]; ++j)
      {

        auto type_j = recv_data[i][mpinf.type * N + j];

        auto pos_x = recv_data[i][(mpinf.pos * N) + (3 * j) + 0];
        auto pos_y = recv_data[i][(mpinf.pos * N) + (3 * j) + 1];
        auto pos_z = recv_data[i][(mpinf.pos * N) + (3 * j) + 2];

        ofs_xyz << type_j;

        if (output_id)
        {
          auto id_j = recv_data[i][mpinf.id * N + j];

          ofs_xyz << " " << id_j;
        }

        ofs_xyz << " " << pos_x + p_o.x << " " << pos_y + p_o.y << " " << pos_z + p_o.z;

        if (output_velocity)
        {
          auto vel_x = recv_data[i][(mpinf.vel * N) + (3 * j) + 0];
          auto vel_y = recv_data[i][(mpinf.vel * N) + (3 * j) + 1];
          auto vel_z = recv_data[i][(mpinf.vel * N) + (3 * j) + 2];
          ofs_xyz << " " << vel_x << " " << vel_y << " " << vel_z;
        }

        if (output_acceleration)
        {
          auto acc_x = recv_data[i][(mpinf.acc * N) + (3 * j) + 0];
          auto acc_y = recv_data[i][(mpinf.acc * N) + (3 * j) + 1];
          auto acc_z = recv_data[i][(mpinf.acc * N) + (3 * j) + 2];
          ofs_xyz << " " << acc_x << " " << acc_y << " " << acc_z;
        }
        ofs_xyz << "\n";
      }
    }

    ofs_xyz << std::flush;
#endif
  }

  //================================================
  //                                              ||
  //================================================
  void Atom_data::dump_xyz_mpi_shared_atoms(int64_t)
  {
#if defined(CAVIAR_WITH_MPI)

    const auto &pos = atom_data->atom_struct_owned.position;
    const auto &vel = atom_data->atom_struct_owned.velocity;
    const auto &acc = atom_data->atom_struct_owned.acceleration;
    const auto &id = atom_data->atom_struct_owned.id;
    const auto &type = atom_data->atom_struct_owned.type;

    unsigned int pos_size = pos.size(); // num_local_atoms;
    const unsigned nprocs = comm->nprocs;
    unsigned int num_total_atoms = pos.size();

    unsigned int num_active_atoms = 0;
    for (unsigned int i = 0; i < pos_size; ++i)
    {    
      if ( atom_data->atom_struct_owned.mpi_rank[i] == my_mpi_rank) num_active_atoms++;
    }
    unsigned int send_num = num_active_atoms;

    std::vector<unsigned int> recv_num(nprocs, 0);
    //==================================// num_local_atoms send
    if (my_mpi_rank != 0)
    {

      MPI_Send(&send_num, 1, MPI::UNSIGNED, 0, 0, MPI::COMM_WORLD);
    }
    if (my_mpi_rank == 0)
    {
      for (unsigned i = 1; i < nprocs; ++i)
      {
        MPI_Recv(&recv_num[i], 1, MPI::UNSIGNED, i, 0, MPI::COMM_WORLD, MPI_STATUS_IGNORE);

        //num_total_atoms += recv_num[i];
      }
    }
    auto &mpinf = mpi_packet_info_owned;

    if (!mpinf.initialized)
    {
      mpinf.initialized = true;
      mpinf.total = 0;
      if (output_id)
      {
        mpinf.id = mpinf.total;
        mpinf.total += 1;
      }
      mpinf.type = mpinf.total;
      mpinf.total += 1;

      mpinf.pos = mpinf.total;
      mpinf.total += 3;

      if (output_velocity)
      {
        mpinf.vel = mpinf.total;
        mpinf.total += 3;
      }
      if (output_acceleration)
      {
        mpinf.acc = mpinf.total;
        mpinf.total += 3;
      }
    }

    std::vector<double> send_data;              // rank != 0
    std::vector<std::vector<double>> recv_data; // rank 0

    if (my_mpi_rank != 0)
    {

      send_data.resize(mpinf.total * send_num, 0);
    }

    if (my_mpi_rank == 0)
    {
      recv_data.resize(nprocs);
      for (unsigned i = 1; i < nprocs; ++i)
      {
        if (recv_num[i] > 0)
          recv_data[i].resize(mpinf.total * recv_num[i], 0);
      }
    }

    if (my_mpi_rank != 0)
    {
      auto N = send_num;

      for (unsigned int i = 0; i < send_num; ++i)
      {
        if ( atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank) continue;

        if (output_id)
          send_data[mpinf.id * N + i] = id[i];

        send_data[mpinf.type * N + i] = type[i];

        send_data[(mpinf.pos * N) + (3 * i) + 0] = pos[i].x;
        send_data[(mpinf.pos * N) + (3 * i) + 1] = pos[i].y;
        send_data[(mpinf.pos * N) + (3 * i) + 2] = pos[i].z;

        if (output_velocity)
        {
          send_data[(mpinf.vel * N) + (3 * i) + 0] = vel[i].x;
          send_data[(mpinf.vel * N) + (3 * i) + 1] = vel[i].y;
          send_data[(mpinf.vel * N) + (3 * i) + 2] = vel[i].z;
        }

        if (output_acceleration)
        {
          send_data[(mpinf.acc * N) + (3 * i) + 0] = acc[i].x;
          send_data[(mpinf.acc * N) + (3 * i) + 1] = acc[i].y;
          send_data[(mpinf.acc * N) + (3 * i) + 2] = acc[i].z;
        }
      }
    }
    //-----------------------------------------------//
    if (my_mpi_rank != 0)
      if (send_num > 0)

        MPI_Send(send_data.data(), mpinf.total * send_num, MPI_DOUBLE, 0, 0, MPI::COMM_WORLD); // TAG 1

    if (my_mpi_rank == 0)
      for (unsigned i = 1; i < nprocs; ++i)
        if (recv_num[i] > 0)

          MPI_Recv(recv_data[i].data(), mpinf.total * recv_num[i], MPI_DOUBLE, i, 0, MPI::COMM_WORLD, MPI_STATUS_IGNORE); // TAG 1

    //-----------------------------------------------//
    if (my_mpi_rank != 0)
      return;

    Vector<double> p_o{0, 0, 0};

    if (position_offset != nullptr)
      p_o = position_offset->current_value;

    ofs_xyz << num_total_atoms << "\nAtom\n";

    for (unsigned int i = 0; i < pos_size; ++i)
    {

      if ( atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank) continue;

      ofs_xyz << type[i];
      if (output_id)
      {

        ofs_xyz << " " << id[i];
      }

      ofs_xyz << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;

      if (output_velocity)
        ofs_xyz << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
      if (output_acceleration)
        ofs_xyz << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
      ofs_xyz << "\n";
    }

    for (unsigned i = 1; i < nprocs; ++i)
    {

      int N = recv_num[i];

      for (unsigned int j = 0; j < recv_num[i]; ++j)
      {

        auto type_j = recv_data[i][mpinf.type * N + j];

        auto pos_x = recv_data[i][(mpinf.pos * N) + (3 * j) + 0];
        auto pos_y = recv_data[i][(mpinf.pos * N) + (3 * j) + 1];
        auto pos_z = recv_data[i][(mpinf.pos * N) + (3 * j) + 2];

        ofs_xyz << type_j;

        if (output_id)
        {
          auto id_j = recv_data[i][mpinf.id * N + j];

          ofs_xyz << " " << id_j;
        }

        ofs_xyz << " " << pos_x + p_o.x << " " << pos_y + p_o.y << " " << pos_z + p_o.z;

        if (output_velocity)
        {
          auto vel_x = recv_data[i][(mpinf.vel * N) + (3 * j) + 0];
          auto vel_y = recv_data[i][(mpinf.vel * N) + (3 * j) + 1];
          auto vel_z = recv_data[i][(mpinf.vel * N) + (3 * j) + 2];
          ofs_xyz << " " << vel_x << " " << vel_y << " " << vel_z;
        }

        if (output_acceleration)
        {
          auto acc_x = recv_data[i][(mpinf.acc * N) + (3 * j) + 0];
          auto acc_y = recv_data[i][(mpinf.acc * N) + (3 * j) + 1];
          auto acc_z = recv_data[i][(mpinf.acc * N) + (3 * j) + 2];
          ofs_xyz << " " << acc_x << " " << acc_y << " " << acc_z;
        }
        ofs_xyz << "\n";
      }
    }

    ofs_xyz << std::flush;
#endif
  }

  //================================================
  //                                              ||
  //================================================

  void Atom_data::dump_xyz_ghost_mpi(int64_t)
  {
#if defined(CAVIAR_WITH_MPI)

    const auto &pos = atom_data->atom_struct_ghost.position;
    const auto &vel = atom_data->atom_struct_ghost.velocity;
    const auto &acc = atom_data->atom_struct_ghost.acceleration;
    const auto &id = atom_data->atom_struct_ghost.id;
    const auto &type = atom_data->atom_struct_ghost.type;

    unsigned int send_num = pos.size(); // num_local_atoms;
    const unsigned nprocs = comm->nprocs;
    unsigned int num_total_atoms = pos.size();

    std::vector<unsigned int> recv_num(nprocs, 0);

    //==================================// num_local_atoms send
    if (my_mpi_rank != 0)
      MPI_Send(&send_num, 1, MPI::UNSIGNED, 0, 0, MPI::COMM_WORLD);

    if (my_mpi_rank == 0)
    {
      for (unsigned i = 1; i < nprocs; ++i)
        MPI_Recv(&recv_num[i], 1, MPI::UNSIGNED, i, 0, MPI::COMM_WORLD, MPI_STATUS_IGNORE);

      for (unsigned int i = 1; i < nprocs; ++i)
        num_total_atoms += recv_num[i];
    }

    auto &mpinf = mpi_packet_info_ghost;

    if (!mpinf.initialized)
    {
      mpinf.initialized = true;
      mpinf.total = 0;

        if (output_id)
    {  
      mpinf.id = mpinf.total;
      mpinf.total += 1;
      }
      mpinf.type = mpinf.total;
      mpinf.total += 1;

      mpinf.pos = mpinf.total;
      mpinf.total += 3;

      if (output_velocity)
      {
        mpinf.vel = mpinf.total;
        mpinf.total += 3;
      }
      if (output_acceleration)
      {
        mpinf.acc = mpinf.total;
        mpinf.total += 3;
      }
    }
    std::vector<double> send_data;              // rank != 0
    std::vector<std::vector<double>> recv_data; // rank 0

    if (my_mpi_rank != 0)
      send_data.resize(mpinf.total * send_num, 0);

    if (my_mpi_rank == 0)
    {
      recv_data.resize(nprocs);
      for (unsigned i = 1; i < nprocs; ++i)
        recv_data[i].resize(mpinf.total * recv_num[i], 0);
    }

    if (my_mpi_rank != 0)
    {
      auto N = send_num;

      for (unsigned int i = 0; i < send_num; ++i)
      {
        if (output_id)
          send_data[mpinf.id * N + i] = id[i];

        send_data[mpinf.type * N + i] = type[i];

        send_data[(mpinf.pos * N) + (3 * i) + 0] = pos[i].x;
        send_data[(mpinf.pos * N) + (3 * i) + 1] = pos[i].y;
        send_data[(mpinf.pos * N) + (3 * i) + 2] = pos[i].z;

        if (output_velocity)
        {
          send_data[(mpinf.vel * N) + (3 * i) + 0] = vel[i].x;
          send_data[(mpinf.vel * N) + (3 * i) + 1] = vel[i].y;
          send_data[(mpinf.vel * N) + (3 * i) + 2] = vel[i].z;
        }

        if (output_acceleration)
        {
          send_data[(mpinf.acc * N) + (3 * i) + 0] = acc[i].x;
          send_data[(mpinf.acc * N) + (3 * i) + 1] = acc[i].y;
          send_data[(mpinf.acc * N) + (3 * i) + 2] = acc[i].z;
        }
      }
    }
    //-----------------------------------------------//

    if (my_mpi_rank != 0)
      if (send_num > 0)
        MPI_Send(send_data.data(), mpinf.total * send_num, MPI_DOUBLE, 0, 0, MPI::COMM_WORLD); // TAG 1

    if (my_mpi_rank == 0)
      for (unsigned i = 1; i < nprocs; ++i)
        if (recv_num[i] > 0)
          MPI_Recv(recv_data[i].data(), mpinf.total * recv_num[i], MPI_DOUBLE, i, 0, MPI::COMM_WORLD, MPI_STATUS_IGNORE); // TAG 1

    //-----------------------------------------------//

    if (my_mpi_rank != 0)
      return;

    Vector<double> p_o{0, 0, 0};

    if (position_offset != nullptr)
      p_o = position_offset->current_value;

    ofs_xyz_ghost << num_total_atoms << "\nAtom\n";

    for (unsigned int i = 0; i < send_num; ++i)
    {
      ofs_xyz_ghost << type[i] ;
      if (output_id)
            ofs_xyz_ghost << " " << id[i] ;

      ofs_xyz_ghost << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;
      if (output_velocity)
        ofs_xyz_ghost << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
      if (output_acceleration)
        ofs_xyz_ghost << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
      ofs_xyz_ghost << "\n";
    }

    for (unsigned i = 1; i < nprocs; ++i)
    {
      int N = recv_num[i];
      for (unsigned int j = 0; j < recv_num[i]; ++j)
      {


        auto type_j = recv_data[i][mpinf.type * N + j];

        auto pos_x = recv_data[i][(mpinf.pos * N) + (3 * j) + 0];
        auto pos_y = recv_data[i][(mpinf.pos * N) + (3 * j) + 1];
        auto pos_z = recv_data[i][(mpinf.pos * N) + (3 * j) + 2];
        ofs_xyz_ghost << type_j ;
        if (output_id)
        {
          auto id_j = recv_data[i][mpinf.id * N + j];

        ofs_xyz_ghost << " " << id_j ;
}
        ofs_xyz_ghost << " " << pos_x + p_o.x << " " << pos_y + p_o.y << " " << pos_z + p_o.z;

        if (output_velocity)
        {
          auto vel_x = recv_data[i][(mpinf.vel * N) + (3 * j) + 0];
          auto vel_y = recv_data[i][(mpinf.vel * N) + (3 * j) + 1];
          auto vel_z = recv_data[i][(mpinf.vel * N) + (3 * j) + 2];
          ofs_xyz_ghost << " " << vel_x << " " << vel_y << " " << vel_z;
        }

        if (output_acceleration)
        {
          auto acc_x = recv_data[i][(mpinf.acc * N) + (3 * j) + 0];
          auto acc_y = recv_data[i][(mpinf.acc * N) + (3 * j) + 1];
          auto acc_z = recv_data[i][(mpinf.acc * N) + (3 * j) + 2];
          ofs_xyz_ghost << " " << acc_x << " " << acc_y << " " << acc_z;
        }
        ofs_xyz_ghost << "\n";
      }
    }

    ofs_xyz_ghost << std::flush;

#endif
  }

  //================================================
  //                                              ||
  //================================================
  void Atom_data::dump_xyz_serial(int64_t, bool mpi_files)
  {
    auto &pos = atom_data->atom_struct_owned.position;
    auto &type = atom_data->atom_struct_owned.type;
    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &acc = atom_data->atom_struct_owned.acceleration;
    auto &id = atom_data->atom_struct_owned.id;

    auto nta = atom_data->atom_struct_owned.position.size();

    Vector<double> p_o{0, 0, 0};
    if (position_offset != nullptr)
      p_o = position_offset->current_value;

    if (!mpi_files)
    {
      ofs_xyz << nta << "\nAtom\n";

      for (unsigned int i = 0; i < nta; ++i)
      {
        ofs_xyz << type[i];
        if (output_id)
          ofs_xyz << " " << id[i];
        ofs_xyz << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;
        if (output_velocity)
          ofs_xyz << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
        if (output_acceleration)
          ofs_xyz << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
        ofs_xyz << "\n";
      }

      ofs_xyz << std::flush;
    }
    else
    {
      unsigned int num_active_atoms = 0;
      for (unsigned int i = 0; i < nta; ++i)
      {    
        if ( atom_data->atom_struct_owned.mpi_rank[i] == my_mpi_rank) num_active_atoms++;
      }

      ofs_xyz_mpi << num_active_atoms << "\nAtom\n";

      for (unsigned int i = 0; i < nta; ++i)
      {
        if ( atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank) continue;

        ofs_xyz_mpi << type[i];
        if (output_id)
          ofs_xyz_mpi << " " << id[i];

        ofs_xyz_mpi << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;
        if (output_velocity)
          ofs_xyz_mpi << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
        if (output_acceleration)
          ofs_xyz_mpi << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
        ofs_xyz_mpi << "\n";
      }

      ofs_xyz_mpi << std::flush;
    }
  }

  //================================================
  //                                              ||
  //================================================
  void Atom_data::dump_xyz_ghost_serial(int64_t, bool mpi_files)
  {
    auto &pos = atom_data->atom_struct_ghost.position;
    auto &type = atom_data->atom_struct_ghost.type;
    auto &vel = atom_data->atom_struct_ghost.velocity;
    auto &acc = atom_data->atom_struct_ghost.acceleration;
    auto &id = atom_data->atom_struct_ghost.id;
    auto nta = atom_data->atom_struct_ghost.position.size();

    Vector<double> p_o{0, 0, 0};
    if (position_offset != nullptr)
      p_o = position_offset->current_value;

    if (!mpi_files)
    {
      ofs_xyz_ghost << nta << "\nAtom\n";

      for (unsigned int i = 0; i < nta; ++i)
      {
        ofs_xyz_ghost << type[i];
        if (output_id)
          ofs_xyz_ghost << " " << id[i];
        ofs_xyz_ghost << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;
        if (output_velocity)
          ofs_xyz_ghost << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
        if (output_acceleration)
          ofs_xyz_ghost << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
        ofs_xyz_ghost << "\n";
      }

      ofs_xyz_ghost << std::flush;
    }
    else
    {
      ofs_xyz_ghost_mpi << nta << "\nAtom\n";

      for (unsigned int i = 0; i < nta; ++i)
      {
        ofs_xyz_ghost_mpi << type[i];
        if (output_id)
          ofs_xyz_ghost_mpi << " " << id[i];

        ofs_xyz_ghost_mpi << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;
        if (output_velocity)
          ofs_xyz_ghost_mpi << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
        if (output_acceleration)
          ofs_xyz_ghost_mpi << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
        ofs_xyz_ghost_mpi << "\n";
      }

      ofs_xyz_ghost_mpi << std::flush;
    }
  }

} // writer

CAVIAR_NAMESPACE_CLOSE
