
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



bool Atom_data::exchange_owned () {
  if (domain==nullptr) error->all("Atom_data::exchange_owned: domain = nullptr");
  bool make_neighlist = false;

  const auto bc = domain->boundary_condition;
  
  // cutoff extra may make problem for induced_charge mesh.
  // in the situations that the mesh size is exactly equal to the domain 
  // measurements, cutoff_extra may makes the particles go outside of the mesh.
  // that throws an exception when calculating the force or accelerations on
  // the particle. in order to avoid that, one can fix it by setting cutoff_extra
  // to zero, or by fixing particle position in the forcefield calculations.
  const auto x_llow = domain->lower_local.x - cutoff_extra;
  const auto x_lupp = domain->upper_local.x + cutoff_extra;
  const auto y_llow = domain->lower_local.y - cutoff_extra;
  const auto y_lupp = domain->upper_local.y + cutoff_extra;
  const auto z_llow = domain->lower_local.z - cutoff_extra;
  const auto z_lupp = domain->upper_local.z + cutoff_extra;

  const auto x_width = domain->upper_local.x - domain->lower_local.x;
  const auto y_width = domain->upper_local.y - domain->lower_local.y;
  const auto z_width = domain->upper_local.z - domain->lower_local.z;


  auto &pos = owned.position;

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  const auto me = domain->me;
  if (me==0) {
  for (unsigned int i = 0 ; i < num_local_atoms; ++i) {
    if (bc.x == 1) { //TODO While should be changed to if in the future
      while (pos[i].x < x_llow) {pos[i].x += x_width; owned.msd_domain_cross[i].x -= 1 ;}
      while (pos[i].x > x_lupp) {pos[i].x -= x_width; owned.msd_domain_cross[i].x += 1 ;}
    }
    if (bc.y == 1) {
      while (pos[i].y < y_llow) {pos[i].y += y_width; owned.msd_domain_cross[i].y -= 1 ;}
      while (pos[i].y > y_lupp) {pos[i].y -= y_width; owned.msd_domain_cross[i].y += 1 ;}      
    }
    if (bc.z == 1) {
      while (pos[i].z < z_llow) {pos[i].z += z_width; owned.msd_domain_cross[i].z -= 1 ;}
      while (pos[i].z > z_lupp) {pos[i].z -= z_width; owned.msd_domain_cross[i].z += 1 ;}      
    }
  } 
  }
#elif defined(CAVIAR_WITH_MPI)
  auto &vel = owned.velocity;
  auto &acc = owned.acceleration;
  auto &id  = owned.id;
  auto &type = owned.type;

  const auto grid_index_x = domain->grid_index_x;
  const auto grid_index_y = domain->grid_index_y;
  const auto grid_index_z = domain->grid_index_z;

  const auto nprocs_x = domain->nprocs_x;
  const auto nprocs_y = domain->nprocs_y;
  const auto nprocs_z = domain->nprocs_z;


  const auto me = domain->me;

  const auto &all = domain->all;

  std::vector<std::vector<std::vector<std::vector<int>>>> o_send_id, o_recv_id, o_send_index;

  o_send_id.resize(3); o_recv_id.resize(3); o_send_index.resize(3);
  for (auto i=0; i<3;++i) {
    o_send_id[i].resize(3);  o_recv_id[i].resize(3); o_send_index[i].resize(3);
    for (auto j=0; j<3;++j) {
      o_send_id[i][j].resize(3); o_recv_id[i][j].resize(3); o_send_index[i][j].resize(3);

    }
  }

  int o_recv_n[3][3][3]; // num of owned to be recieved from domain [i][j][k]

  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        o_recv_n[i][j][k] = 0;
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

    int x_val = 0, y_val = 0, z_val = 0;

    if (xlc) x_val = -1;
    if (xuc) x_val = +1;
    if (ylc) y_val = -1;
    if (yuc) y_val = +1;
    if (zlc) z_val = -1;
    if (zuc) z_val = +1;


    if (grid_index_x == 0) x_val *= bc.x;             // periodic or non-periodic boundary condition
    if (grid_index_x == nprocs_x - 1) x_val *= bc.x;  // //
    if (grid_index_y == 0) y_val *= bc.y;              // //
    if (grid_index_y == nprocs_y - 1) y_val *= bc.y;  // //
    if (grid_index_z == 0) z_val *= bc.z;              // //
    if (grid_index_z == nprocs_z - 1) z_val *= bc.z;  // //

    if (x_val == 0 && y_val == 0 && z_val == 0) {
      continue;
    } else { 
      o_send_id    [x_val + 1][y_val + 1][z_val + 1].emplace_back (id[i]);
      o_send_index [x_val + 1][y_val + 1][z_val + 1].emplace_back (i);
    }

  }

  MPI_Barrier (mpi_comm);

// ===============================send num of owned


  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me != all[i][j][k]) {
          int num_s = o_send_id[i][j][k].size();
          MPI_Send (&num_s, 1, MPI_INT, all[i][j][k], 0, mpi_comm);
          MPI_Recv (&o_recv_n[i][j][k], 1, MPI_INT, all[i][j][k], 0, mpi_comm, MPI_STATUS_IGNORE);
        }
      }
    }
  }

  MPI_Barrier (mpi_comm);

// ===============================send id of owned
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me != all[i][j][k]) {
          int num_s = o_send_id[i][j][k].size();  
          if (num_s>0){
            MPI_Send (o_send_id[i][j][k].data(), num_s, MPI_INT, all[i][j][k] , 1, mpi_comm);
          }          
          unsigned num_r = o_recv_n [i][j][k];
          if (num_r>0){
            int *tmp_r = new int[num_r]; // TODO maybe change it to a static vector
            MPI_Recv (&tmp_r, num_r, MPI_INT, all[i][j][k], 1, mpi_comm, MPI_STATUS_IGNORE);
            for (unsigned m = 0; m < num_r; ++m)
              o_recv_id[i][j][k].push_back (tmp_r[m]);
            delete[] tmp_r;
          }
        }
      }
    }
  }

  MPI_Barrier (mpi_comm);
// ==============================send type of owned
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me != all[i][j][k]) {
          for (auto m : o_send_index[i][j][k]) {
            MPI_Send (&type[m], 1, MPI_INT, all[i][j][k], id[m], mpi_comm);
          }
          for (auto m : o_recv_id[i][j][k]) {
            int tmp;
            MPI_Recv (&tmp, 1, MPI_INT, all[i][j][k], m, mpi_comm, MPI_STATUS_IGNORE);
            type.emplace_back (tmp);
            id.emplace_back (m);
          } 
        }
      }
    }
  }

  MPI_Barrier (mpi_comm);
// ==============================send position of owned
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me != all[i][j][k]) {
          for (auto m : o_send_index[i][j][k]) { 
            Vector<double> p_tmp {pos[m].x, pos[m].y, pos[m].z};
            MPI_Send (&p_tmp.x, 3, MPI_DOUBLE, all[i][j][k], id[m], mpi_comm);
          }
          for (auto m : o_recv_id[i][j][k]) {
            Vector<double> p_tmp;
            MPI_Recv (&p_tmp, 3, MPI_DOUBLE, all[i][j][k], m, mpi_comm, MPI_STATUS_IGNORE);
            if (i==0) while (p_tmp.x < x_llow) p_tmp.x+= x_width;
            if (j==0) while (p_tmp.y < y_llow) p_tmp.y+= y_width;
            if (k==0) while (p_tmp.z < z_llow) p_tmp.z+= z_width;
            if (i==2) while (p_tmp.x > x_lupp) p_tmp.x-= x_width;
            if (j==2) while (p_tmp.y > y_lupp) p_tmp.y-= y_width;
            if (k==2) while (p_tmp.z > z_lupp) p_tmp.z-= z_width;
            pos.emplace_back (p_tmp.x, p_tmp.y, p_tmp.z);
            ++num_local_atoms;
          } 
        }
      }
    }
  }
  MPI_Barrier (mpi_comm);
// ==============================send velocity of owned
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me != all[i][j][k]) {
          for (auto m : o_send_index[i][j][k]) { 
            Vector<double> p_tmp {vel[m].x, vel[m].y, vel[m].z}; 
            MPI_Send (&p_tmp.x, 3, MPI_DOUBLE, all[i][j][k], id[m], mpi_comm); 
          } 
          for (auto m : o_recv_id[i][j][k]) { 
            Vector<double> p_tmp; 
            MPI_Recv (&p_tmp, 3, MPI_DOUBLE, all[i][j][k], m, mpi_comm, MPI_STATUS_IGNORE); 
            vel.emplace_back (p_tmp.x, p_tmp.y, p_tmp.z);
          }  
        }
      }
    }
  }
  MPI_Barrier (mpi_comm);
// ==============================send acceleration of owned
  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
        if (me != all[i][j][k]) {
          for (auto m : o_send_index[i][j][k]) { 
            Vector<double> p_tmp {acc[m].x, acc[m].y, acc[m].z}; 
            MPI_Send (&p_tmp.x, 3, MPI_DOUBLE, all[i][j][k], id[m], mpi_comm); 
          } 
          for (auto m : o_recv_id[i][j][k]) { 
            Vector<double> p_tmp; 
            MPI_Recv (&p_tmp, 3, MPI_DOUBLE, all[i][j][k], m, mpi_comm, MPI_STATUS_IGNORE); 
            acc.emplace_back (p_tmp.x, p_tmp.y, p_tmp.z);
          }  
        }
      }
    }
  }
  MPI_Barrier (mpi_comm);

// ================================================ self move

  std::vector<int>o_send_index_lin; // gather all the indexes in a one dimensional array. used in erase.

  for (auto i=0; i<3;++i) {
    for (auto j=0; j<3;++j) {
      for (auto k=0; k<3;++k) {
//        if (i==1&&j==1&&k==1) continue;
        if (me == all[i][j][k]) {
          auto ii=i-1,jj=j-1,kk=k-1;
          for (auto m : o_send_index[i][j][k]) {
            pos[m].x -= ii * x_width;
            pos[m].y -= jj * y_width;
            pos[m].z -= kk * z_width;            
          }  
        } else {
          for (auto m : o_send_index[i][j][k]) {
            o_send_index_lin.emplace_back (m);
          }
        }
      }
    }
  }
  MPI_Barrier (mpi_comm);


// ================================================ delete moved atoms
  if (o_send_index_lin.size()>0) {
    remove_atom(o_send_index_lin);
    make_neighlist = true;   
  }


#else

#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for
#endif   
  for (unsigned int i = 0 ; i < num_local_atoms; ++i) {
    if (bc.x == 1) { 
      while (pos[i].x < x_llow) {pos[i].x += x_width; owned.msd_domain_cross[i].x -= 1 ;}
      while (pos[i].x > x_lupp) {pos[i].x -= x_width; owned.msd_domain_cross[i].x += 1 ;}
    }
    if (bc.y == 1) {
      while (pos[i].y < y_llow) {pos[i].y += y_width; owned.msd_domain_cross[i].y -= 1 ;}
      while (pos[i].y > y_lupp) {pos[i].y -= y_width; owned.msd_domain_cross[i].y += 1 ;}      
    }
    if (bc.z == 1) {
      while (pos[i].z < z_llow) {pos[i].z += z_width; owned.msd_domain_cross[i].z -= 1 ;}
      while (pos[i].z > z_lupp) {pos[i].z -= z_width; owned.msd_domain_cross[i].z += 1 ;}      
    }
  }  
#endif
  return make_neighlist;
}



CAVIAR_NAMESPACE_CLOSE


