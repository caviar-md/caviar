
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

#include "caviar/objects/md_simulator.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/integrator.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/force_field.h"
#include "caviar/objects/constraint.h"
#include "caviar/objects/writer.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/interpreter/communicator.h"

namespace caviar {

namespace objects {

Md_simulator::Md_simulator (CAVIAR *fptr) : Pointers{fptr},
    atom_data{nullptr}, integrator{nullptr},
    initialized{false} {
  FC_OBJECT_INITIALIZE
}

Md_simulator::~Md_simulator () {}

void Md_simulator::verify_settings () {
  
}


bool Md_simulator::read_base_class_commands (caviar::interpreter::Parser *) {
  bool command_called = false;
/*
  parser -> keep_current_token();
  auto token = parser -> get_val_token();
  auto t = token.string_value;
  if (string_cmp(t,"output_total_force") ) {
    command_called = true;
    output_total_force();
  } else if (string_cmp(t,"output_total_energy") ) {
    command_called = true;
    output_total_energy();
  }
*/
  return command_called;
}

void Md_simulator::initialize () {
  initialized = true;
}

void Md_simulator::step (int i) {

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


    integrator -> step_part_I ();


    for (auto&& c : constraint)
      c -> step_part_I (i);


    integrator -> step_part_II ();



    for (auto&& c : constraint)
      c -> step_part_II (i);

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

#else


  integrator -> step_part_I ();


  for (auto&& c : constraint)
    c -> step_part_I (i);


  integrator -> step_part_II ();

  for (auto&& c : constraint)
      c -> step_part_II (i);

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


  integrator->step_part_III ();

  for (auto&& c : constraint)
    c -> step_part_III (i);

  for (auto&& w : writer) 
    w -> write (i, time);

#endif

  time += dt;

}


bool Md_simulator::boundary_condition () {
  bool result = false;

  // the result only in MPI case can
  // become true, due to exchanging a particle between domains
  result = atom_data -> exchange_owned (); 

  atom_data -> exchange_ghost ();
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#elif defined(CAVIAR_WITH_MPI)
  MPI_Allreduce (MPI::IN_PLACE, &result, 1, MPI::BOOL, MPI::LOR, mpi_comm);
#endif
  return result;
}

void Md_simulator::setup () {
  t_start = clock();
  if (use_time) {
    initial_step = (int) initial_time/dt;
    final_step = (int) final_time/dt;
  }

  FC_NULLPTR_CHECK(atom_data)
  atom_data->record_owned_old_data();

  if (neighborlist.size()==0) output->warning("Md_simulator::setup: neighborlist.size() = 0");
  if (force_field.size()==0) output->warning("Md_simulator::setup: force_field.size() = 0");

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



void Md_simulator::cleanup () {

}


}

} // namespace caviar

