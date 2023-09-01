
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

  const auto bc = domain->boundary_condition;

  // cutoff extra may make problem for induced_charge mesh.
  // in the situations that the mesh size is exactly equal to the domain
  // measurements, cutoff_extra may makes the particles go outside of the mesh.
  // that throws an exception when calculating the force or accelerations on
  // the particle. in order to avoid that, one can fix it by setting cutoff_extra
  // to zero, or by fixing particle position in the forcefield calculations.
  // const auto x_llow = domain->lower_global.x - cutoff_extra;
  // const auto x_lupp = domain->upper_global.x + cutoff_extra;
  // const auto y_llow = domain->lower_global.y - cutoff_extra;
  // const auto y_lupp = domain->upper_global.y + cutoff_extra;
  // const auto z_llow = domain->lower_global.z - cutoff_extra;
  // const auto z_lupp = domain->upper_global.z + cutoff_extra;

  const auto x_llow = domain->lower_global.x;
  const auto x_lupp = domain->upper_global.x;
  const auto y_llow = domain->lower_global.y;
  const auto y_lupp = domain->upper_global.y;
  const auto z_llow = domain->lower_global.z;
  const auto z_lupp = domain->upper_global.z;

  const auto x_width = domain->upper_global.x - domain->lower_global.x;
  const auto y_width = domain->upper_global.y - domain->lower_global.y;
  const auto z_width = domain->upper_global.z - domain->lower_global.z;

  auto &pos = atom_struct_owned.position;

  unsigned int num_local_atoms = atom_struct_owned.id.size();

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < num_local_atoms; ++i)
  {
    if (bc.x == 1)
    { 
      if (pos[i].x < x_llow) // while or if
      {
        pos[i].x += x_width;
        if (msd_process)
          atom_struct_owned.msd_domain_cross[i].x -= 1;
      }
      else if (pos[i].x > x_lupp)
      {
        pos[i].x -= x_width;
        if (msd_process)
          atom_struct_owned.msd_domain_cross[i].x += 1;
      }
    }
    if (bc.y == 1)
    {
      if (pos[i].y < y_llow)
      {
        pos[i].y += y_width;
        if (msd_process)

          atom_struct_owned.msd_domain_cross[i].y -= 1;
      }
      else if (pos[i].y > y_lupp)
      {
        pos[i].y -= y_width;
        if (msd_process)

          atom_struct_owned.msd_domain_cross[i].y += 1;
      }
    }
    if (bc.z == 1)
    {
      if (pos[i].z < z_llow)
      {
        pos[i].z += z_width;
        if (msd_process)

          atom_struct_owned.msd_domain_cross[i].z -= 1;
      }
      else if (pos[i].z > z_lupp)
      {
        pos[i].z -= z_width;
        if (msd_process)

          atom_struct_owned.msd_domain_cross[i].z += 1;
      }
    }
  }
  return update_verlet_list;
}
//======================================================
//                                                    ||
//                                                    ||
//======================================================

CAVIAR_NAMESPACE_CLOSE
