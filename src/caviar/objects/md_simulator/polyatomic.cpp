
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

#include "caviar/objects/md_simulator/polyatomic.h"
#include "caviar/objects/integrator.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/force_field.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/constraint.h"
#include "caviar/objects/writer.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/interpreter/communicator.h"
#include <ctime>

namespace caviar {
namespace objects {
namespace md_simulator {

Polyatomic::Polyatomic (CAVIAR *fptr) : Md_simulator{fptr}
{
  FC_OBJECT_INITIALIZE_INFO
  use_time = false;
  use_step = false;
}

Polyatomic::~Polyatomic () {}

bool Polyatomic::read (caviar::interpreter::Parser *parser) {
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

bool Polyatomic::run () {
  output->info ("MD started.");
  verify_settings ();

  t_start = clock();

  for (auto i = initial_step; i < final_step; ++i) 
    step(i);

  cleanup(); 

  output->info ("MD finished.");

  t_end = clock();
  double simulation_time =  ( (float)t_end - (float)t_start ) / CLOCKS_PER_SEC;


#if defined (CAVIAR_WITH_MPI) 
  std::string s = "process " + std::to_string(comm->me) + " : ";
#else
  std::string s = "";
#endif
  s += "simulation time: " + std::to_string(simulation_time) + " (seconds)";
  std::cout << s << std::endl;
  //output->info (s, 3);
  return true; //WARNING
}

void Polyatomic::verify_settings () {
  FC_NULLPTR_CHECK(atom_data)
  FC_NULLPTR_CHECK(integrator)
  if (use_time && use_step)
    error->all ("simulator md: verify_settings: cannot use both 'time' and 'step'.");
  if (!(use_time || use_step))
    error->all ("simulator md: verify_settings: please assign 'time' or 'step' variables.");
}




void Polyatomic::step (int i) {

  FC_OBJECT_VERIFY_SETTINGS

  if (!initialized) {
    initialize();

    for (auto&& c : constraint)
      c -> verify_settings ();

    time = dt*initial_step;
    setup();
  }

  atom_data->record_owned_old_data();

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  static int my_mpi_rank = comm->me;

  if (my_mpi_rank==0) {
   
    bool make_list = boundary_condition (); 


    if (make_list) 
      for (auto&& n: neighborlist) n->build_neighlist ();
    else {
      for (auto&& n: neighborlist)
        make_list = make_list || n->rebuild_neighlist ();
      if (make_list) 
        for (auto&& n: neighborlist) n->build_neighlist ();
    }

  }


  atom_data -> synch_owned_data(0);

  atom_data -> reset_owned_acceleration();

  for (auto&& f : force_field)
    f -> calculate_acceleration ();

  MPI_Barrier(mpi_comm);

  if (my_mpi_rank==0) {


    integrator -> step_part_III ();


    for (auto&& c : constraint)
      c -> step_part_III (i);
  }


  for (auto&& w : writer) 
    w -> write (i, time);


  if (my_mpi_rank==0) {


    integrator -> step_part_I ();


    for (auto&& c : constraint)
      c -> step_part_I (i);


    integrator -> step_part_II ();



    for (auto&& c : constraint)
      c -> step_part_II (i);
  }


#else


  bool make_list = boundary_condition (); 


  if (make_list) 
    for (auto&& n: neighborlist) n->build_neighlist ();
  else {
    for (auto&& n: neighborlist)
      make_list = make_list || n->rebuild_neighlist ();
    if (make_list) 
      for (auto&& n: neighborlist) n->build_neighlist ();
  }

  atom_data -> reset_owned_acceleration();

  for (auto&& f : force_field)
    f -> calculate_acceleration ();


  integrator -> step_part_I ();


  for (auto&& c : constraint)
    c -> step_part_I (i);


  integrator -> step_part_II ();

  for (auto&& c : constraint)
      c -> step_part_II (i);


  integrator->step_part_III ();

  for (auto&& c : constraint)
    c -> step_part_III (i);


  for (auto&& w : writer) 
    w -> write (i, time);
#endif

  time += dt;

}


/*
void Polyatomic::setup () {
  t_start = clock();
  if (use_time) {
    initial_step = (int) initial_time/dt;
    final_step = (int) final_time/dt;
  }

  FC_NULLPTR_CHECK(atom_data)
  atom_data->record_owned_old_data();

  if (neighborlist.size()==0) output->warning("Polyatomic::setup: neighborlist.size() = 0");
  if (force_field.size()==0) output->warning("Polyatomic::setup: force_field.size() = 0");

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

  atom_data->synch_owned_data(0);

  //if (my_mpi_rank==0) { // XXX STRANGE MPI BUG! if uncommented, the code will stuck somewhere else.

    for (auto&& n: neighborlist) n->init ();

    atom_data -> exchange_owned (); // isn't neccesary

    atom_data -> exchange_ghost ();

    for (auto&& n: neighborlist)
      n->build_neighlist ();

  //}


  atom_data -> reset_owned_acceleration();

  for (auto&& f : force_field) 
    f -> calculate_acceleration ();


  for (auto&& w : writer) 
    w -> write (0, 0.0);

#else

  for (auto&& n: neighborlist) n->init ();

  atom_data -> exchange_owned (); // isn't neccesary

  atom_data -> exchange_ghost ();

  for (auto&& n: neighborlist) n->build_neighlist ();

  atom_data -> reset_owned_acceleration();

  for (auto&& f : force_field) 
    f -> calculate_acceleration ();

  for (auto&& w : writer) 
    w -> write (0, 0.0);

#endif

}
*/



} //md_simulator
} //objects
} // namespace caviar

