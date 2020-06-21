
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
#include "caviar/utility/python_utils_def.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"
#include "caviar/utility/interpreter_io_headers.h"
#include <ctime>
#include <sys/stat.h> // used for mkdir()

namespace caviar {
namespace objects {
namespace writer {

Atom_data::Atom_data (CAVIAR *fptr) : Writer{fptr},
    atom_data{nullptr}, domain{nullptr},
    energy_step{100}, xyz_step{1000}, povray_step{10000}, msd_step{1000},
    output_energy{false}, output_xyz{false}, output_povray{false},
    output_msd{false},
    output_velocity{false},  output_acceleration{false}
{
  FC_OBJECT_INITIALIZE_INFO
  msd_initial_step = 0;
  tStart1 = clock();
}

Atom_data::~Atom_data () {
  if (ofs_xyz.is_open())    ofs_xyz.close();
  if (ofs_energy.is_open()) ofs_energy.close();
  if (ofs_povray.is_open()) ofs_povray.close();
  if (ofs_msd.is_open())    ofs_msd.close();
}

bool Atom_data::read (caviar::interpreter::Parser *) {
  /*
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"set_domain") || string_cmp(t,"domain")) {
      FIND_OBJECT_BY_NAME(domain,it)
      domain = object_container->domain[it->second.index];
    } else  if (string_cmp(t,"xyz_step")) {
      GET_OR_CHOOSE_A_INT(xyz_step,"","")
      output_xyz = true;
    } else  if (string_cmp(t,"povray_step")) {
      GET_OR_CHOOSE_A_INT(povray_step,"","")
      output_povray = true;
    } else  if (string_cmp(t,"energy_step")) {
      GET_OR_CHOOSE_A_INT(energy_step,"","")
      output_energy = true;
    } else  if (string_cmp(t,"msd_step")) {
      GET_OR_CHOOSE_A_INT(msd_step,"","")
      output_msd = true;
    } else  if (string_cmp(t,"msd_initial_step")) {
      GET_OR_CHOOSE_A_INT(msd_initial_step,"","")
      output_msd = true;
    } else  if (string_cmp(t,"open_files")) {
      open_files();
    } else  if (string_cmp(t,"close_files")) {
      close_files();
    } else  if (string_cmp(t,"output_velocity")) {
      output_velocity = true;
    } else  if (string_cmp(t,"output_acceleration")) {
      //output_acceleration = true;
      std::ofstream ofs ("o_acc");
      const auto &pos = atom_data -> owned.position;  
      const auto &acc = atom_data -> owned.acceleration;  
      for (unsigned int i=0;i<pos.size();++i) {
        ofs << i << " " << acc[i].x << "\t" << acc[i].y << "\t" << acc[i].z << "\n" ;
      }
    } 
    else FC_ERR_UNDEFINED_VAR(t)
  }
  return in_file;
  */
  return true;
}

void Atom_data::initialize(){
  initialized = true;

  FC_NULLPTR_CHECK(atom_data)

  if(output_msd) {
    FC_NULLPTR_CHECK(domain)
  }

// --- just to make povray outpuy folder ---
  if (my_mpi_rank==0) {
    if (output_povray) {
      std::string str_folder_pov;
      str_folder_pov.append("o_pov"); 
      const char* char_folder_pov = str_folder_pov.c_str();
      mkdir (char_folder_pov,0777);  // make povray  output folder //
    }
  }
/*
#ifdef CAVIAR_WITH_MPI
  sprintf ( buffer, "_me%u", comm->me );
  str_filename.append ( buffer);
#endif
*/

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
if (my_mpi_rank==0) {
  if (output_xyz) 
    if (!ofs_xyz.is_open()) 
      ofs_xyz.open("o_xyz.xyz");

  if (output_povray)     
    if (!ofs_povray.is_open()) 
      ofs_povray.open("o_pov.pov");

  if (output_energy)     
    if (!ofs_energy.is_open()) 
      ofs_energy.open("o_energy.txt");

  if (output_msd)     
    if (!ofs_msd.is_open()) 
      ofs_msd.open("o_msd.txt");
}
#elif defined(CAVIAR_WITH_MPI)
  if (my_mpi_rank==0) {
  if (output_xyz) 
    if (!ofs_xyz.is_open()) 
      ofs_xyz.open("o_xyz.xyz");

  if (output_povray)     
    if (!ofs_povray.is_open()) 
      ofs_povray.open("o_pov.pov");

  if (output_energy)     
    if (!ofs_energy.is_open()) 
      ofs_energy.open("o_energy.txt");

  if (output_msd)     
    if (!ofs_msd.is_open()) 
      ofs_msd.open("o_msd.txt");
  }//*/
#else

  if (output_xyz) 
    if (!ofs_xyz.is_open()) {
      ofs_xyz.open("o_xyz.xyz");

    }

  if (output_povray)     
    if (!ofs_povray.is_open()) 
      ofs_povray.open("o_pov.pov");

  if (output_energy)     
    if (!ofs_energy.is_open()) 
      ofs_energy.open("o_energy.txt");

  if (output_msd)     
    if (!ofs_msd.is_open()) 
      ofs_msd.open("o_msd.txt");
#endif
  
}

void Atom_data::open_files(){}
void Atom_data::close_files(){}
void Atom_data::generate(){}
void Atom_data::write(){}
void Atom_data::write(int){} // current time_step
void Atom_data::write(double){} // current time

void Atom_data::write(int i, double t){

  if (!initialized) 
    initialize();

  if (output_xyz) 
    if (i%xyz_step==0)
      dump_xyz(i);
  
  if (output_energy) 
    if (i%energy_step==0)
      dump_energy(i);

  if (output_povray) 
    if (i%povray_step==0)
      dump_povray(i);

  if (output_msd) 
    if (i%msd_step==0)
      dump_msd(i, t);
}

void Atom_data::start_new_files(){} //add_time_to_previous
void Atom_data::start_new_files(std::string &){} //add_time_to_previous



FC_PYDEF_SETGET_PTR(Atom_data,atom_data,objects::Atom_data);
FC_PYDEF_SETGET_PTR(Atom_data,domain,objects::Domain);

/*
FC_PYDEF_SETGET_STDVEC2D(Atom_data,epsilon,Real_t);  
FC_PYDEF_SETGET_STDVEC2D(Atom_data,sigma,Real_t);
FC_PYDEF_SETGET_STDVEC(Atom_data,epsilon_atom,Real_t);  
FC_PYDEF_SETGET_STDVEC(Atom_data,sigma_atom,Real_t);
FC_PYDEF_SETGET_STDVEC2D(Atom_data,cutoff_list,Real_t);
*/
void export_py_Atom_data () {

  using namespace boost::python;

  implicitly_convertible<std::shared_ptr<writer::Atom_data>,          
                         std::shared_ptr<Writer> >(); 




  class_<writer::Atom_data,boost::noncopyable>("Atom_data",init<caviar::CAVIAR*>())

    .def("open_files",&writer::Atom_data::open_files)      
    .def("close_files",&writer::Atom_data::close_files)      


    .def_readwrite("povray_step",&writer::Atom_data::povray_step)      
    .def_readwrite("output_povray",&writer::Atom_data::output_povray)            

    .def_readwrite("xyz_step",&writer::Atom_data::xyz_step)
    .def_readwrite("output_xyz",&writer::Atom_data::output_xyz)            

    .def_readwrite("energy_step",&writer::Atom_data::energy_step)            
    .def_readwrite("output_energy",&writer::Atom_data::output_energy)            

    .def_readwrite("msd_step",&writer::Atom_data::msd_step)            
    .def_readwrite("msd_initial_step",&writer::Atom_data::msd_initial_step)            
    .def_readwrite("output_msd",&writer::Atom_data::output_msd)            

    .def_readwrite("output_velocity",&writer::Atom_data::output_velocity)            

    .add_property("atom_data", &writer::Atom_data::get_atom_data, &writer::Atom_data::set_atom_data)
    .add_property("domain", &writer::Atom_data::get_domain, &writer::Atom_data::set_domain)


  ;
}


/* //XXX
    } else  if (string_cmp(t,"output_acceleration")) {
      //output_acceleration = true;
      std::ofstream ofs ("o_acc");
      const auto &pos = atom_data -> owned.position;  
      const auto &acc = atom_data -> owned.acceleration;  
      for (unsigned int i=0;i<pos.size();++i) {
        ofs << i << " " << acc[i].x << "\t" << acc[i].y << "\t" << acc[i].z << "\n" ;
      }
    } 
*/

} //Atom_data
} //objects
} // namespace caviar

