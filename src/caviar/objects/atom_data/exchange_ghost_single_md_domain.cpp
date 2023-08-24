
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

void Atom_data::exchange_ghost_single_md_domain(long) // timestep
{
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  if (domain->me != 0)
    return;
#endif

  atom_struct_ghost.position.clear();
  atom_struct_ghost.velocity.clear();
  atom_struct_ghost.id.clear();
  atom_struct_ghost.type.clear();

  const auto bc = domain->boundary_condition;

  const auto x_llow = domain->lower_local.x + ghost_cutoff;
  const auto x_lupp = domain->upper_local.x - ghost_cutoff;
  const auto y_llow = domain->lower_local.y + ghost_cutoff;
  const auto y_lupp = domain->upper_local.y - ghost_cutoff;
  const auto z_llow = domain->lower_local.z + ghost_cutoff;
  const auto z_lupp = domain->upper_local.z - ghost_cutoff;

  const auto x_width = domain->upper_local.x - domain->lower_local.x;
  const auto y_width = domain->upper_local.y - domain->lower_local.y;
  const auto z_width = domain->upper_local.z - domain->lower_local.z;

  auto &pos = atom_struct_owned.position;
  auto &vel = atom_struct_owned.velocity;
  auto &id = atom_struct_owned.id;
  auto &type = atom_struct_owned.type;

  auto &g_pos = atom_struct_ghost.position;
  auto &g_vel = atom_struct_ghost.velocity;
  auto &g_id = atom_struct_ghost.id;
  auto &g_type = atom_struct_ghost.type;

  unsigned int num_local_atoms = atom_struct_owned.id.size();

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < num_local_atoms; ++i)
  {
    const auto xlc = pos[i].x < x_llow;
    const auto xuc = pos[i].x > x_lupp;
    const auto ylc = pos[i].y < y_llow;
    const auto yuc = pos[i].y > y_lupp;
    const auto zlc = pos[i].z < z_llow;
    const auto zuc = pos[i].z > z_lupp;

    int x_val, y_val, z_val;
    if (xlc)
      x_val = -1;
    else if (xuc)
      x_val = +1;
    else
      x_val = 0;

    if (ylc)
      y_val = -1;
    else if (yuc)
      y_val = +1;
    else
      y_val = 0;

    if (zlc)
      z_val = -1;
    else if (zuc)
      z_val = +1;
    else
      z_val = 0;

    x_val *= bc.x;
    y_val *= bc.y;
    z_val *= bc.z; // boundary condition

    // not sure if this 'make_ghost_velocity' condition makes much change in the
    // serial code or for low number of particles.
    if (make_ghost_velocity)
    {
      if (x_val != 0)
      {
        g_pos.emplace_back(pos[i].x - x_val * x_width, pos[i].y, pos[i].z);
        g_vel.emplace_back(vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (y_val != 0)
      {
        g_pos.emplace_back(pos[i].x, pos[i].y - y_val * y_width, pos[i].z);
        g_vel.emplace_back(vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (z_val != 0)
      {
        g_pos.emplace_back(pos[i].x, pos[i].y, pos[i].z - z_val * z_width);
        g_vel.emplace_back(vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (x_val != 0 && y_val != 0)
      {
        g_pos.emplace_back(pos[i].x - x_val * x_width, pos[i].y - y_val * y_width, pos[i].z);
        g_vel.emplace_back(vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (x_val != 0 && z_val != 0)
      {
        g_pos.emplace_back(pos[i].x - x_val * x_width, pos[i].y, pos[i].z - z_val * z_width);
        g_vel.emplace_back(vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (y_val != 0 && z_val != 0)
      {
        g_pos.emplace_back(pos[i].x, pos[i].y - y_val * y_width, pos[i].z - z_val * z_width);
        g_vel.emplace_back(vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (x_val != 0 && y_val != 0 && z_val != 0)
      {
        g_pos.emplace_back(pos[i].x - x_val * x_width, pos[i].y - y_val * y_width, pos[i].z - z_val * z_width);
        g_vel.emplace_back(vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
    }
    else
    {
      if (x_val != 0)
      {
        g_pos.emplace_back(pos[i].x - x_val * x_width, pos[i].y, pos[i].z);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (y_val != 0)
      {
        g_pos.emplace_back(pos[i].x, pos[i].y - y_val * y_width, pos[i].z);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (z_val != 0)
      {
        g_pos.emplace_back(pos[i].x, pos[i].y, pos[i].z - z_val * z_width);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (x_val != 0 && y_val != 0)
      {
        g_pos.emplace_back(pos[i].x - x_val * x_width, pos[i].y - y_val * y_width, pos[i].z);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (x_val != 0 && z_val != 0)
      {
        g_pos.emplace_back(pos[i].x - x_val * x_width, pos[i].y, pos[i].z - z_val * z_width);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (y_val != 0 && z_val != 0)
      {
        g_pos.emplace_back(pos[i].x, pos[i].y - y_val * y_width, pos[i].z - z_val * z_width);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
      if (x_val != 0 && y_val != 0 && z_val != 0)
      {
        g_pos.emplace_back(pos[i].x - x_val * x_width, pos[i].y - y_val * y_width, pos[i].z - z_val * z_width);
        g_id.emplace_back(id[i]);
        g_type.emplace_back(type[i]);
      }
    }
  }
}


CAVIAR_NAMESPACE_CLOSE
