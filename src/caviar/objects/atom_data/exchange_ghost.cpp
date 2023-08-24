
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

CAVIAR_NAMESPACE_OPEN


void Atom_data::exchange_ghost(long i) // timestep
{
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

  exchange_ghost_single_md_domain();

#elif defined(CAVIAR_WITH_MPI)

  exchange_ghost_mpi_shared_atoms(i);
  // exchange_ghost_mpi(i);

#else

  exchange_ghost_single_md_domain(i);

#endif
}

//  if (self_ghost_check())
//    error->all (FC_FILE_LINE_FUNC_PARSE, "Self ghost can happen. Force field cutoff is larger than half of a domain.");
/*
bool ::self_ghost_check () {
  const auto x_llow = domain->lower_local.x;
  const auto x_lupp = domain->upper_local.x;
  const auto y_llow = domain->lower_local.y;
  const auto y_lupp = domain->upper_local.y;
  const auto z_llow = domain->lower_local.z;
  const auto z_lupp = domain->upper_local.z;

  const auto x_width = x_lupp - x_llow;
  const auto y_width = y_lupp - y_llow;
  const auto z_width = z_lupp - z_llow;

  const auto cutoff = force_field->cutoff;
  if (2*cutoff>x_width || 2*cutoff>y_width || 2*cutoff>z_width)
    return true;
  return false;
}*/

CAVIAR_NAMESPACE_CLOSE
