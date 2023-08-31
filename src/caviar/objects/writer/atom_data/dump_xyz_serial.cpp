
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
