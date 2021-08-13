
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

namespace caviar {

namespace objects {

void Atom_data::exchange_ghost () {

  ghost.position.clear ();
  ghost.velocity.clear ();
  ghost.id.clear ();
  ghost.type.clear ();
  ghost_rank.clear ();

  
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


  auto &pos = owned.position;
  auto &vel = owned.velocity;
  auto &id  = owned.id;
  auto &type = owned.type;

  auto &g_pos = ghost.position;
  auto &g_vel = ghost.velocity;
  auto &g_id  = ghost.id;
  auto &g_type = ghost.type;


#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  const auto me = domain->me;
  if (me==0) {
  for (unsigned int i=0; i<num_local_atoms; ++i) {
    const auto xlc = pos[i].x < x_llow;
    const auto xuc = pos[i].x > x_lupp;
    const auto ylc = pos[i].y < y_llow;
    const auto yuc = pos[i].y > y_lupp;
    const auto zlc = pos[i].z < z_llow;
    const auto zuc = pos[i].z > z_lupp;

    int x_val, y_val, z_val;
    if (xlc) x_val = -1;
    else if (xuc) x_val = +1;
    else x_val = 0;

    if (ylc) y_val = -1;
    else if (yuc) y_val = +1;
    else y_val = 0;

    if (zlc) z_val = -1;
    else if (zuc) z_val = +1;
    else z_val = 0;

  
    x_val *= bc.x; y_val *= bc.y; z_val *= bc.z; // boundary condition

    // not sure if this 'make_ghost_velocity' condition makes much change in the
    // serial code or for low number of particles.
    if (make_ghost_velocity) {
    if (x_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y, pos[i].z);
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (y_val != 0) { 
        g_pos.emplace_back (pos[i].x, pos[i].y - y_val*y_width, pos[i].z);
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (z_val != 0) {
        g_pos.emplace_back (pos[i].x, pos[i].y, pos[i].z - z_val*z_width);
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && y_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y - y_val*y_width, pos[i].z);
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && z_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y, pos[i].z - z_val*z_width);
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (y_val != 0 && z_val != 0) {
        g_pos.emplace_back (pos[i].x, pos[i].y - y_val*y_width, pos[i].z - z_val*z_width);
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && y_val != 0 && z_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y - y_val*y_width, pos[i].z - z_val*z_width);
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    } else {
    if (x_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y, pos[i].z);
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (y_val != 0) { 
        g_pos.emplace_back (pos[i].x, pos[i].y - y_val*y_width, pos[i].z);
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (z_val != 0) {
        g_pos.emplace_back (pos[i].x, pos[i].y, pos[i].z - z_val*z_width);
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && y_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y - y_val*y_width, pos[i].z);
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && z_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y, pos[i].z - z_val*z_width);
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (y_val != 0 && z_val != 0) {
        g_pos.emplace_back (pos[i].x, pos[i].y - y_val*y_width, pos[i].z - z_val*z_width);
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && y_val != 0 && z_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y - y_val*y_width, pos[i].z - z_val*z_width);
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    }
  }
  }
#elif defined(CAVIAR_WITH_MPI)
  const auto grid_index_x = domain->grid_index_x;
  const auto grid_index_y = domain->grid_index_y;
  const auto grid_index_z = domain->grid_index_z;

  const auto nprocs_x = domain->nprocs_x;
  const auto nprocs_y = domain->nprocs_y;
  const auto nprocs_z = domain->nprocs_z;

//  const auto nprocs = domain->nprocs;
  const auto me = domain->me;
  const auto &all = domain->all;

  std::vector<std::vector<std::vector<std::vector<int>>>> g_send_id, g_recv_id, g_send_index;

  g_send_id.resize(3); g_recv_id.resize(3); g_send_index.resize(3);
  for (auto i=0; i<3;++i) {
    g_send_id[i].resize(3);  g_recv_id[i].resize(3); g_send_index[i].resize(3);
    for (auto j=0; j<3;++j) {
      g_send_id[i][j].resize(3); g_recv_id[i][j].resize(3); g_send_index[i][j].resize(3);

    }
  }


  int g_recv_n[3][3][3]; // num of ghosts to be recieved from domain [i][j][k]

  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        g_recv_n[i][j][k] = 0;
      }
    }
  }

  for (unsigned i=0; i<num_local_atoms; ++i) {
    const auto xlc = pos[i].x < x_llow;
    const auto xuc = pos[i].x > x_lupp;
    const auto ylc = pos[i].y < y_llow;
    const auto yuc = pos[i].y > y_lupp;
    const auto zlc = pos[i].z < z_llow;
    const auto zuc = pos[i].z > z_lupp;

    int x_val, y_val, z_val;
    if (xlc) x_val = -1;
    else if (xuc) x_val = +1;
    else x_val = 0;

    if (ylc) y_val = -1;
    else if (yuc) y_val = +1;
    else y_val = 0;

    if (zlc) z_val = -1;
    else if (zuc) z_val = +1;
    else z_val = 0;
    
    if (grid_index_x == 0) x_val *= bc.x;             // periodic or non-periodic boundary condition
    if (grid_index_x == nprocs_x - 1) x_val *= bc.x;  // //
    if (grid_index_y == 0) y_val *= bc.y;              // //
    if (grid_index_y == nprocs_y - 1) y_val *= bc.y;  // //
    if (grid_index_z == 0) z_val *= bc.z;              // //
    if (grid_index_z == nprocs_z - 1) z_val *= bc.z;  // //

    if (x_val != 0) { 
      g_send_id    [x_val + 1][1][1].emplace_back (id[i]);
      g_send_index [x_val + 1][1][1].emplace_back (i);
    }
    if (y_val != 0) { 
      g_send_id    [1][y_val + 1][1].emplace_back (id[i]);
      g_send_index [1][y_val + 1][1].emplace_back (i);
    }
    if (z_val != 0) { 
      g_send_id    [1][1][z_val + 1].emplace_back (id[i]);
      g_send_index [1][1][z_val + 1].emplace_back (i);
    }
    if (x_val != 0 && y_val != 0) { 
      g_send_id    [x_val + 1][y_val + 1][1].emplace_back (id[i]);
      g_send_index [x_val + 1][y_val + 1][1].emplace_back (i);
    }
    if (x_val != 0 && z_val != 0) { 
      g_send_id    [x_val + 1][1][z_val + 1].emplace_back (id[i]);
      g_send_index [x_val + 1][1][z_val + 1].emplace_back (i);
    }
    if (y_val != 0 && z_val != 0) { 
      g_send_id    [1] [y_val + 1][z_val + 1].emplace_back (id[i]);
      g_send_index [1] [y_val + 1][z_val + 1].emplace_back (i);
    }
    if (x_val != 0 && y_val != 0 && z_val != 0) { 
      g_send_id    [x_val + 1][y_val + 1][z_val + 1].emplace_back (id[i]);
      g_send_index [x_val + 1][y_val + 1][z_val + 1].emplace_back (i);
    }
  }
  MPI_Barrier (mpi_comm);
// ===============================send num of ghosts
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me != all[i][j][k]) {
          int num_s = g_send_id[i][j][k].size();
          MPI_Send (&num_s, 1, MPI_INT, all[i][j][k], 0, mpi_comm);
          MPI_Recv (&g_recv_n[i][j][k], 1, MPI_INT, all[i][j][k], 0, mpi_comm, MPI_STATUS_IGNORE);
        }
      }
    }
  }
  MPI_Barrier (mpi_comm);
// ===============================send id of ghosts
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me != all[i][j][k]) {
          int num_s = g_send_id[i][j][k].size();  
          if (num_s>0){
            MPI_Send (g_send_id[i][j][k].data(), num_s, MPI_INT, all[i][j][k] , 1, mpi_comm);
          }          
          unsigned num_r = g_recv_n [i][j][k];
          if (num_r>0){
            int *tmp_r = new int[num_r];
            MPI_Recv (tmp_r, num_r, MPI_INT, all[i][j][k], 1, mpi_comm, MPI_STATUS_IGNORE);
            for (unsigned m = 0; m < num_r; ++m)
              g_recv_id[i][j][k].emplace_back (tmp_r[m]);
            delete[] tmp_r;
          }
        }
      }
    }
  }
  MPI_Barrier (mpi_comm);
// ==============================send type of  ghosts
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me != all[i][j][k]) {
          for (auto m : g_send_index[i][j][k]) {
            MPI_Send (&type[m], 1, MPI_INT, all[i][j][k], id[m], mpi_comm);
          }
          for (auto m : g_recv_id[i][j][k]) {
            int tmp;
            MPI_Recv (&tmp, 1, MPI_INT, all[i][j][k], m, mpi_comm, MPI_STATUS_IGNORE);
            g_type.emplace_back (tmp);
            g_id.emplace_back (m);
            ghost_rank.emplace_back (all[i][j][k]);
          } 
        }
      }
    }
  }
  MPI_Barrier (mpi_comm);
// ==============================send position of ghosts
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me != all[i][j][k]) {          
          for (auto m : g_send_index[i][j][k]) { 
            Vector<double> p_tmp {pos[m].x, pos[m].y, pos[m].z};
            MPI_Send (&p_tmp.x, 3, MPI_DOUBLE, all[i][j][k], id[m], mpi_comm);
          }
          for (auto m : g_recv_id[i][j][k]) {
            Vector<double> p_tmp;
            MPI_Recv (&p_tmp, 3, MPI_DOUBLE, all[i][j][k], m, mpi_comm, MPI_STATUS_IGNORE);
            if (i==0) while (p_tmp.x < x_llow) p_tmp.x+= x_width;
            if (j==0) while (p_tmp.y < y_llow) p_tmp.y+= y_width;
            if (k==0) while (p_tmp.z < z_llow) p_tmp.z+= z_width;
            if (i==2) while (p_tmp.x > x_lupp) p_tmp.x-= x_width;
            if (j==2) while (p_tmp.y > y_lupp) p_tmp.y-= y_width;
            if (k==2) while (p_tmp.z > z_lupp) p_tmp.z-= z_width;
            g_pos.emplace_back (p_tmp.x, p_tmp.y, p_tmp.z);
          } 
        }
      }
    }
  }
  MPI_Barrier (mpi_comm);
// ==============================send velocity of ghosts
  if (make_ghost_velocity) {
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
          if (me != all[i][j][k]) {
          for (auto m : g_send_index[i][j][k]) { 
            Vector<double> p_tmp {vel[m].x, vel[m].y, vel[m].z}; 
            MPI_Send (&p_tmp.x, 3, MPI_DOUBLE, all[i][j][k], id[m], mpi_comm); 
          } 
          for (auto m : g_recv_id[i][j][k]) { 
            Vector<double> p_tmp; 
            MPI_Recv (&p_tmp, 3, MPI_DOUBLE, all[i][j][k], m, mpi_comm, MPI_STATUS_IGNORE);
            g_vel.emplace_back (p_tmp.x, p_tmp.y, p_tmp.z); 
          }  
        }
      }
    }
  }
  MPI_Barrier (mpi_comm);
  }
  //===========================SAME DOMAIN GHOSTS
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me == all[i][j][k]) {
          auto x_val=i-1,y_val=j-1,z_val=k-1;
          for (auto m :g_send_index[i][j][k]) {
            g_vel.emplace_back (vel[m].x, vel[m].y, vel[m].z);
            g_pos.emplace_back (pos[m].x - x_val*x_width, pos[m].y - y_val*y_width, pos[m].z - z_val*z_width); 
            g_id.emplace_back (id[m]); 
            g_type.emplace_back (type[m]);
            ghost_rank.emplace_back (me);
          }    
        }
      }
    }
  }
#else

  for (unsigned int i=0; i<num_local_atoms; ++i) {
    const auto xlc = pos[i].x < x_llow;
    const auto xuc = pos[i].x > x_lupp;
    const auto ylc = pos[i].y < y_llow;
    const auto yuc = pos[i].y > y_lupp;
    const auto zlc = pos[i].z < z_llow;
    const auto zuc = pos[i].z > z_lupp;

    int x_val, y_val, z_val;
    if (xlc) x_val = -1;
    else if (xuc) x_val = +1;
    else x_val = 0;

    if (ylc) y_val = -1;
    else if (yuc) y_val = +1;
    else y_val = 0;

    if (zlc) z_val = -1;
    else if (zuc) z_val = +1;
    else z_val = 0;

  
    x_val *= bc.x; y_val *= bc.y; z_val *= bc.z; // boundary condition

    // not sure if this 'make_ghost_velocity' condition makes much change in the
    // serial code or for low number of particles.
    if (make_ghost_velocity) {
    if (x_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y, pos[i].z); 
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (y_val != 0) { 
        g_pos.emplace_back (pos[i].x, pos[i].y - y_val*y_width, pos[i].z); 
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (z_val != 0) {
        g_pos.emplace_back (pos[i].x, pos[i].y, pos[i].z - z_val*z_width);
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && y_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y - y_val*y_width, pos[i].z);
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && z_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y, pos[i].z - z_val*z_width); 
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (y_val != 0 && z_val != 0) {
        g_pos.emplace_back (pos[i].x, pos[i].y - y_val*y_width, pos[i].z - z_val*z_width); 
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && y_val != 0 && z_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y - y_val*y_width, pos[i].z - z_val*z_width); 
        g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]); 
        g_type.emplace_back (type[i]);
    }
    } else {
    if (x_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y, pos[i].z); 
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (y_val != 0) { 
        g_pos.emplace_back (pos[i].x, pos[i].y - y_val*y_width, pos[i].z); 
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (z_val != 0) {
        g_pos.emplace_back (pos[i].x, pos[i].y, pos[i].z - z_val*z_width);
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && y_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y - y_val*y_width, pos[i].z);
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && z_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y, pos[i].z - z_val*z_width); 
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (y_val != 0 && z_val != 0) {
        g_pos.emplace_back (pos[i].x, pos[i].y - y_val*y_width, pos[i].z - z_val*z_width); 
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]);
        g_type.emplace_back (type[i]);
    }
    if (x_val != 0 && y_val != 0 && z_val != 0) { 
        g_pos.emplace_back (pos[i].x - x_val*x_width, pos[i].y - y_val*y_width, pos[i].z - z_val*z_width); 
        //g_vel.emplace_back (vel[i].x, vel[i].y, vel[i].z);
        g_id.emplace_back (id[i]); 
        g_type.emplace_back (type[i]);
    }
    }
  }

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





} //objects

} // namespace caviar


