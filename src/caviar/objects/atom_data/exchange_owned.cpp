
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


CAVIAR_NAMESPACE_OPEN


bool Atom_data::exchange_owned(long i) // timestep
{
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

  return exchange_owned_single_md_domain();

#elif defined(CAVIAR_WITH_MPI)

  return exchange_owned_mpi_shared_atoms(i);
  // return exchange_owned_mpi(i);

#else

  return exchange_owned_single_md_domain(i);

#endif
}

CAVIAR_NAMESPACE_CLOSE
