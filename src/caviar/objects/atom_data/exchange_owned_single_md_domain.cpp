
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

//======================================================
//                                                    ||
//                                                    ||
//======================================================

bool Atom_data::exchange_owned_single_md_domain(long) // timestep
{
  if (domain == nullptr)
    error->all("Atom_data::exchange_owned: domain = nullptr");

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  if (domain->me != 0)
    return false;
#endif

  bool update_verlet_list = false;

  auto &pos = atom_struct_owned.position;

  auto pos_size = pos.size();

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < pos_size; ++i)
  {
    caviar::Vector<int> msd {0,0,0};
    pos[i] = domain->fix_position(pos[i], msd, update_verlet_list);

    if (msd_process)
      atom_struct_owned.msd_domain_cross[i] +=  msd;
  }
  return update_verlet_list;
}
//======================================================
//                                                    ||
//                                                    ||
//======================================================

CAVIAR_NAMESPACE_CLOSE
