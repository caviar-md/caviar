
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
#include "caviar/interpreter/communicator.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace writer
{
  void Atom_data::dump_povray_mpi(int64_t, double)
  {

  }

  void Atom_data::dump_povray_mpi_shared_atoms(int64_t, double)
  {
    
  }

  void Atom_data::dump_povray_serial(int64_t, double)
  {
    /*
    #ifdef CAVIAR_WITH_MPI
      MPI_Barrier (mpi_comm);

      unsigned nprocs = comm->nprocs;
      const int me = comm->me;

      const auto &pos = atom_struct_owned.position;
      const auto &id  = atom_struct_owned.id;
      const auto &type = atom_struct_owned.type;

      std::vector<std::vector<int>> all_id;
      std::vector<Vector<double>> all_pos;
      std::vector<int> all_type;

      const auto nla = pos.size();//atom_data -> num_local_atoms;
      const auto nta =  num_total_atoms;
    //  std::cout <<"nla" << nla <<std::endl;
    //  std::cout <<"nta" << nta <<std::endl;
      unsigned *nla_list = new unsigned [ nprocs ];

      if (nprocs > 1) {


        for (auto j=0;j<atom_data->atom_struct_ghost.id.size(); ++j)
          if (atom_data->atom_struct_ghost.id[j] == 26)
            std::cout<<"me:" << me << "\tt:" << i <<std::endl;

        for (auto j=0;j<pos.size(); ++j)
          if (id[j] == 26)
            std::cout<<"me:" << me << " p:" << atom_data->atom_struct_owned.acceleration[j] << "\tt:" << i << "\ttyp:" << type[j] <<std::endl;
        std::cout<<"me:" << me << " g:" << atom_data->atom_struct_ghost.position.size() << std::endl;
        std::cout<<"me:" << me << " " << pos.size() << id.size() << type.size() << atom_data->atom_struct_owned.acceleration.size() << atom_data -> atom_struct_owned.velocity.size() <<std::endl;

        if (me != 0) { //==================================// nla send
          MPI_Send (&nla, 1, MPI::UNSIGNED, 0, 0, mpi_comm);
        } else {
          nla_list[0]=nla;
          for (unsigned int i=1; i < nprocs; ++i) {
            MPI_Recv (nla_list+i, 1, MPI::UNSIGNED, i, 0, mpi_comm, MPI_STATUS_IGNORE);
          }
        } //-----------------------------------------------//

        MPI_Barrier (mpi_comm);

        if (me != 0) { //==================================// id send
          MPI_Send (id.data(), nla, MPI::UNSIGNED, 0, 0, mpi_comm); //  CHANGE id to INT before send it  //
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


        if (me != 0) { //==================================// pos send
          for (unsigned int i=0; i<pos.size (); ++i){
            Vector<double> p_tmp {pos[i].x, pos[i].y, pos[i].z};
            MPI_Send (&p_tmp.x, 3, MPI_DOUBLE, 0, id[i], mpi_comm);
          }
        } else {
          all_pos.reserve(nta);
          for (unsigned int i=0; i<nla; ++i){
            Vector<double> p_tmp;
            p_tmp = pos[i];
            all_pos.push_back (p_tmp);
          }

          for (unsigned int i=1; i < nprocs; ++i) {
            for (auto j : all_id[i]) {
              Vector<double> p_tmp;
              MPI_Recv (&p_tmp, 3, MPI_DOUBLE, i, j, mpi_comm, MPI_STATUS_IGNORE);
              all_pos.push_back (p_tmp);
            }
          }
        } //-----------------------------------------------//

        MPI_Barrier (mpi_comm);

        if (me != 0) { //==================================// type send
          for (unsigned int i=0; i<pos.size (); ++i){
            MPI_Send (&type[i], 1, MPI_INT, 0, id[i], mpi_comm);
          }
        } else {
          all_type.reserve(nta);
          for (unsigned int i=0; i<nla; ++i){
            int tmp = type[i];
            all_type.push_back (tmp);
          }

          for (unsigned i=1; i < nprocs; ++i) {
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
      delete[] nla_list;
      MPI_Barrier (mpi_comm);
    #endif

    #ifdef CAVIAR_WITH_MPI
      if (me == 0) {
    #endif


      std::string str_povray = "o_pov/o_";

      int b = i;int c = 0;
      while (b>0) {
        b /= 10;
        ++c;
      }

      int d = 10 - c;
      while (d>0) {
        str_povray += std::to_string (0);
        --d;
      }

      str_povray += std::to_string (i);
      str_povray.append (".pov");

      const char * char_povray  = str_povray.c_str ();

      ofs_povray.open (char_povray);


      const auto &all_pos = atom_data -> atom_struct_owned.position;
      //const auto &all_type = atom_data ->atom_struct_owned.type;
      unsigned int nta = atom_data -> num_total_atoms;
      ofs_povray << "\nunion {\n";
      for (unsigned int i = 0; i<nta; ++i) {
        ofs_povray << "\tsphere {"
                 << "<" << all_pos[i].x << "," << all_pos[i].y
                << "," << all_pos[i].z << ">,"
    //            << geometry->radius[all_type[i]-1] << "}\n";
                << "1.0" << "}\n";
      }
      ofs_povray << "\ttexture {\n"
                 << "\t\tpigment { color Yellow }\n }\n";

      ofs_povray <<"}\n";

      ofs_povray.close ();

    #ifdef CAVIAR_WITH_MPI
      }
    #endif
    */
  }

} // writer

CAVIAR_NAMESPACE_CLOSE
