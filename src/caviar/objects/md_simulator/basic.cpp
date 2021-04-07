
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

#include "caviar/objects/md_simulator/basic.h"
#include "caviar/objects/integrator.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/utility/time_utility.h"
#include "caviar/interpreter/communicator.h"
#ifdef CAVIAR_WITH_OPENMP
#include <omp.h>
#endif

#include <ctime>

namespace caviar {
namespace objects {
namespace md_simulator {

Basic::Basic (CAVIAR *fptr) : Md_simulator{fptr}
{
  FC_OBJECT_INITIALIZE_INFO
  use_time = false;
  use_step = false;
}

Basic::~Basic () {}

bool Basic::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"integrator") || string_cmp(t,"set_integrator")) {
      FIND_OBJECT_BY_NAME(integrator,it)
      integrator = object_container->integrator[it->second.index];
    } else if (string_cmp(t,"dt")) {
      GET_OR_CHOOSE_A_REAL(dt,"","")
    } else if (string_cmp(t,"initial_time")) {
      GET_OR_CHOOSE_A_REAL(initial_time,"","")
      use_time = true;
    } else if (string_cmp(t,"final_time")) {
      GET_OR_CHOOSE_A_REAL(final_time,"","")
      use_time = true;
    } else if (string_cmp(t,"initial_step")) {
      GET_OR_CHOOSE_A_INT(initial_step,"","")
      use_step = true;
    } else if (string_cmp(t,"final_step")) {
      GET_OR_CHOOSE_A_INT(final_step,"","")
      use_step = true;
    } else if (string_cmp(t,"run")) {
      run();
    } else if (string_cmp(t,"add_constraint") || string_cmp(t,"constraint")) {
      FIND_OBJECT_BY_NAME(constraint,it)
      constraint.push_back(object_container->constraint[it->second.index]);
    } else if (string_cmp(t,"add_writer") || string_cmp(t,"writer")) {
      FIND_OBJECT_BY_NAME(writer,it)
      writer.push_back(object_container->writer[it->second.index]);
    } else if (string_cmp(t,"add_neighborlist") || string_cmp(t,"neighborlist")) {
      FIND_OBJECT_BY_NAME(neighborlist,it)
      neighborlist.push_back(object_container->neighborlist[it->second.index]);
    } else if (string_cmp(t,"add_force_field") || string_cmp(t,"force_field")) {
      FIND_OBJECT_BY_NAME(force_field,it)
      force_field.push_back(object_container->force_field[it->second.index]);
    } else  if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else  if (string_cmp(t,"step")) {
      int i=0;
      GET_OR_CHOOSE_A_INT(i,"","")
      step (i);
      continue;
    } else if (read_base_class_commands(parser)) {
    } else error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
  }
  return in_file;
}

bool Basic::run () {
  output->info ("MD started.");
#ifdef CAVIAR_WITH_OPENMP
  output->info("CAVIAR is started with OpenMP. Max number of threads is " + std::to_string(omp_get_max_threads()));
#endif  
  verify_settings ();
  
  double wall0 = get_wall_time();
  //double cpu0  = get_cpu_time();
      
  //t_start = clock();

  for (auto i = initial_step; i < final_step; ++i) 
    step(i);

  double wall1 = get_wall_time();
  //double cpu1  = get_cpu_time();
      
  cleanup(); 

  output->info ("MD finished.");

  //t_end = clock();
  //double simulation_time =  ( (float)t_end - (float)t_start ) / CLOCKS_PER_SEC;


  // #if defined (CAVIAR_WITH_MPI) 
  //   std::string s = "MPI process number" + std::to_string(comm->me) + " : ";
  // #else
  //   std::string s = "";
  // #endif  
  
  std::string s = "md_simulation run call:" ;
  s += "Wall Time = " + std::to_string(wall1 - wall0) + " sec. , ";  
  //s += "CPU Time  = " + std::to_string(cpu1  - cpu0) + " sec. ";                      
   
  output->info (s);
  
  return true; //WARNING
}

void Basic::verify_settings () {
  FC_NULLPTR_CHECK(atom_data)
  FC_NULLPTR_CHECK(integrator)
  if (use_time && use_step)
    error->all ("simulator md: verify_settings: cannot use both 'time' and 'step'.");
  if (!(use_time || use_step))
    error->all ("simulator md: verify_settings: please assign 'time' or 'step' variables.");
}

} //md_simulator
} //objects
} // namespace caviar

