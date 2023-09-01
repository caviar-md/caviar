
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
        // dump_xyz_mpi(i);
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
        dump_xyz_ghost_serial(i, true);
    }
    else
    {
      dump_xyz_ghost_serial(i);
    }

#else

    dump_xyz_ghost_serial(i);

#endif
  }

} // writer

CAVIAR_NAMESPACE_CLOSE
