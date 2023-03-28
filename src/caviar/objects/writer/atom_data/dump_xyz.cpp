
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

namespace caviar {


namespace objects {
namespace writer {







void Atom_data::dump_xyz (int64_t i) {

//#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#if defined(CAVIAR_WITH_MPI)
  MPI_Barrier (mpi_comm);  


  const auto &pos = atom_data -> owned.position;
  //const auto &vel = atom_data -> owned.velocity;
  //const auto &acc = atom_data -> owned.acceleration;
  const auto &id  = atom_data -> owned.id; 
  const auto &type = atom_data ->owned.type;  

  std::vector<std::vector<int>> all_id;
  std::vector<Vector<double>> all_pos;
  //std::vector<Vector<double>> all_vel;
  //std::vector<Vector<double>> all_acc;

  std::vector<int> all_type;
  
  
  const auto nla = pos.size();//atom_data -> num_local_atoms;
  const auto nta = atom_data -> num_total_atoms;
  const unsigned nprocs = comm -> nprocs;
  unsigned *nla_list = new unsigned [nprocs];

  if (nprocs > 1) {

    if (my_mpi_rank != 0) { //==================================// nla send
      MPI_Send (&nla, 1, MPI::UNSIGNED, 0, 0, mpi_comm);
    } else {
      nla_list[0]=nla;
      for (unsigned i=1; i < nprocs; ++i) {
        MPI_Recv (nla_list+i, 1, MPI::UNSIGNED, i, 0, mpi_comm, MPI_STATUS_IGNORE);
      }
    } //-----------------------------------------------//

    MPI_Barrier (mpi_comm);

    if (my_mpi_rank != 0) { //==================================// id send
      MPI_Send (id.data(), nla, MPI::UNSIGNED, 0, 0, mpi_comm); // CHANGE id to INT before send it //
    } else {

      all_id.resize (nprocs);

      for (unsigned int i=0; i < nla; ++i) { // local copy
        all_id[0].push_back (id[i]);
      }
      for (unsigned int i=1; i < nprocs; ++i) {
        unsigned *tmp = new unsigned [nla_list[i]];
        MPI_Recv (tmp, nla_list[i], MPI::UNSIGNED, i, 0, mpi_comm, MPI_STATUS_IGNORE); // XXX tmp or &tmp
        for (unsigned int j=0; j<nla_list[i];++j)
          all_id[i].push_back(tmp[j]);
        delete[] tmp;
      }
    } //-----------------------------------------------//

      MPI_Barrier (mpi_comm);


    if (my_mpi_rank != 0) { //==================================// pos send
      for (unsigned i=0; i<pos.size (); ++i){
        Vector<double> p_tmp {pos[i].x, pos[i].y, pos[i].z};
        MPI_Send (&p_tmp.x, 3, MPI_DOUBLE, 0, id[i], mpi_comm);
      }
    } else {
      all_pos.reserve(nta);
      for (unsigned i=0; i<nla; ++i){
        Vector<double> p_tmp;
        p_tmp = pos[i];
        all_pos.push_back (p_tmp);
      }

      for (unsigned i=1; i < nprocs; ++i) {
        for (auto j : all_id[i]) {
          Vector<double> p_tmp;
          MPI_Recv (&p_tmp, 3, MPI_DOUBLE, i, j, mpi_comm, MPI_STATUS_IGNORE);
          all_pos.push_back (p_tmp);
        }
      }
    } //-----------------------------------------------//

    MPI_Barrier (mpi_comm);

    if (my_mpi_rank != 0) { //==================================// type send
      for (unsigned int i=0; i<pos.size (); ++i){
        MPI_Send (&type[i], 1, MPI_INT, 0, id[i], mpi_comm); 
      }
    } else {
      all_type.reserve(nta);
      for (unsigned int i=0; i<nla; ++i){
        int tmp = type[i];
        all_type.push_back (tmp);
      }

      for (unsigned int i=1; i < nprocs; ++i) {
        for (auto j : all_id[i]) { 
          int tmp;
          MPI_Recv (&tmp, 1, MPI_INT, i, j, mpi_comm, MPI_STATUS_IGNORE);
          all_type.push_back (tmp);
        }
      }
    } //-----------------------------------------------//

    MPI_Barrier (mpi_comm);


  } else {//==================================// one_processor MPI counterpart.
          // this part can be written in a more concise way.
      all_id.resize (nprocs);
      for (unsigned int i=0; i < nla; ++i) { 
        all_id[0].push_back (id[i]);
        all_pos.push_back (pos[i]);
        all_type.push_back (type[i]);

      }
  
  } //-----------------------------------------------//

  MPI_Barrier (mpi_comm);

  Vector<double> p_o {0,0,0};
  if (position_offset != nullptr) p_o = position_offset->current_value;

  if (my_mpi_rank == 0) {

    ofs_xyz << all_type.size() << "\nAtom\n";

    for (unsigned int i=0; i<all_type.size(); ++i) {
      ofs_xyz << all_type[i] << " " << all_pos[i].x + p_o.x<< " " << all_pos[i].y + p_o.y<< " " << all_pos[i].z + p_o.z << "\n";
    }


    ofs_xyz << std::flush;
  }

  delete[] nla_list;

#else


  auto &all_pos = atom_data -> owned.position;
  auto &all_type = atom_data ->owned.type;  
  auto &all_vel = atom_data -> owned.velocity;
  auto &all_acc = atom_data -> owned.acceleration;
  auto nta = atom_data -> owned.position.size();

  ofs_xyz << nta << "\nAtom\n";

  //if (my_mpi_rank==0) {

  Vector<double> p_o {0,0,0};
  if (position_offset != nullptr) p_o = position_offset->current_value;

  if (output_velocity && output_acceleration) {
    for (unsigned int i = 0; i<nta; ++i) 
    {
      ofs_xyz << all_type[i] << " " << all_pos[i].x + p_o.x<< " " << all_pos[i].y + p_o.y<< " " << all_pos[i].z + p_o.z ;
      ofs_xyz << " " << all_vel[i].x << " " << all_vel[i].y << " " << all_vel[i].z ;
      ofs_xyz << " " << all_acc[i].x << " " << all_acc[i].y << " " << all_acc[i].z ;
      ofs_xyz << "\n";
    }
  } else if (output_velocity) {
    for (unsigned int i = 0; i<nta; ++i) 
    {
      ofs_xyz << all_type[i] << " " << all_pos[i].x + p_o.x<< " " << all_pos[i].y + p_o.y<< " " << all_pos[i].z + p_o.z ;
      ofs_xyz << " " << all_vel[i].x << " " << all_vel[i].y << " " << all_vel[i].z ;
      ofs_xyz << "\n";
    }
  } else if (output_acceleration) {
    for (unsigned int i = 0; i<nta; ++i) 
    {
      ofs_xyz << all_type[i] << " " << all_pos[i].x + p_o.x<< " " << all_pos[i].y + p_o.y<< " " << all_pos[i].z + p_o.z ;
      ofs_xyz << " " << all_acc[i].x << " " << all_acc[i].y << " " << all_acc[i].z ;
      ofs_xyz << "\n";
    }
  } else {
    for (unsigned int i = 0; i<nta; ++i) 
    {
      ofs_xyz << all_type[i] << " " << all_pos[i].x + p_o.x<< " " << all_pos[i].y + p_o.y<< " " << all_pos[i].z + p_o.z ;
      ofs_xyz << "\n";
    }
  }

  ofs_xyz << std::flush;
  //}
#endif
  
  double wallTimeXyzDump2 = get_wall_time();      
  
  double dtstart= wallTimeXyzDump2 - wallTimeXyzDump1;
  std::string s = "writer::atom_data:: dump_xyz at step " + std::to_string(i) +
                + " . Elapsed time since previous xyz dump: " + std::to_string(dtstart);
  output->info (s, 2);
  wallTimeXyzDump1 = wallTimeXyzDump2;
}


} // writer 
} //objects 

} // namespace caviar

