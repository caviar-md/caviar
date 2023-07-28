
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

void Atom_data::synch_owned_data(int)
{
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

  auto mpi_fc_vector_type = comm->mpi_fc_vector_type;
  auto my_mpi_rank = comm->me;
  if (comm->nprocs == 1)
    return;

  int root_rank = 0;

  int pos_size = atom_struct_owned.position.size();
  int root_pos_size = pos_size;

  int type_size = atom_struct_owned.type.size();
  // int root_type_size = type_size;

  // int mass_size = atom_type_params.mass.size();
  // int root_mass_size = mass_size;

  // int charge_size = atom_type_params.charge.size();
  // int root_charge_size = charge_size;

  MPI_Bcast(&root_pos_size, 1, MPI::INT, 0, mpi_comm);

  // XXX not necessary yet.
  // MPI_Bcast (&synch_owned_data_bcast_details,   1, MPI::BOOL, 0, mpi_comm);

  // if (synch_owned_data_bcast_details)
  // {
  //   MPI_Bcast(&root_type_size, 1, MPI::INT, 0, mpi_comm);
  //   MPI_Bcast(&root_mass_size, 1, MPI::INT, 0, mpi_comm);
  //   MPI_Bcast(&root_charge_size, 1, MPI::INT, 0, mpi_comm);
  // }

  if (my_mpi_rank != root_rank)
  {

    if (pos_size != root_pos_size)
    {
      atom_struct_owned.position.resize(root_pos_size, caviar::Vector<double>{0, 0, 0});
      atom_struct_owned.velocity.resize(root_pos_size, caviar::Vector<double>{0, 0, 0});
      atom_struct_owned.acceleration.resize(root_pos_size, caviar::Vector<double>{0, 0, 0});
    }

    if (type_size != root_pos_size)
    {
      atom_struct_owned.type.resize(root_type_size, 0);
    }

    // if (mass_size != root_mass_size)
    // {
    //   atom_type_params.mass.resize(root_mass_size, 1.0);
    // }

    // if (charge_size != root_charge_size)
    // {
    //   atom_type_params.charge.resize(root_charge_size, 0.0);
    // }
  }

  MPI_Bcast(&atom_struct_owned.position[0], root_pos_size, mpi_fc_vector_type, 0, mpi_comm);

  // if (synch_owned_data_bcast_details)
  // {
  //   MPI_Bcast(&atom_struct_owned.type[0], root_type_size, MPI::INT, 0, mpi_comm);

  //   MPI_Bcast(&atom_type_params.mass[0], root_mass_size, MPI::DOUBLE, 0, mpi_comm);

  //   MPI_Bcast(&atom_type_params.charge[0], root_charge_size, MPI::DOUBLE, 0, mpi_comm);
  // }

  synch_owned_data_bcast_details = false;

#endif
}

CAVIAR_NAMESPACE_CLOSE
