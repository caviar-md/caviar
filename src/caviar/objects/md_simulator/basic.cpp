
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
#include "caviar/utility/python_utils_def.h"
#include "caviar/objects/integrator.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/interpreter/communicator.h"
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

bool Basic::read (caviar::interpreter::Parser *) {
  /*
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
  */
  return true;
}

bool Basic::run () {
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

void Basic::verify_settings () {
  FC_NULLPTR_CHECK(atom_data)
  FC_NULLPTR_CHECK(integrator)
  if (use_time && use_step)
    error->all ("simulator md: verify_settings: cannot use both 'time' and 'step'.");
  if (!(use_time || use_step))
    error->all ("simulator md: verify_settings: please assign 'time' or 'step' variables.");
}


void Basic::add_neighborlist (const std::shared_ptr<objects::Neighborlist > &obj) {
  neighborlist.emplace_back(obj);
}

void Basic::add_force_field (const std::shared_ptr<objects::Force_field > &obj) {
  force_field.emplace_back(obj);
}

void Basic::add_constraint (const std::shared_ptr<objects::Constraint > &obj) {
  constraint.emplace_back(obj);
}

void Basic::add_writer (const std::shared_ptr<objects::Writer > &obj) {
  writer.emplace_back(obj);
}


FC_PYDEF_SETGET_PTR(Basic,atom_data,Atom_data);
FC_PYDEF_SETGET_PTR(Basic,integrator,Integrator);


  double time, dt;
  clock_t t_start, t_end;
  double initial_time, final_time;
  bool use_time, use_step;
  int initial_step, final_step;
  bool initialized;

void export_py_Basic () {

  using namespace boost::python;

  implicitly_convertible<std::shared_ptr<md_simulator::Basic>,          
                         std::shared_ptr<Md_simulator> >(); 

  class_<md_simulator::Basic,boost::noncopyable>("Basic",init<caviar::CAVIAR*>())
    .def("run",&md_simulator::Basic::run)
    .def("add_neighborlist",&md_simulator::Basic::add_neighborlist)
    .def("add_force_field",&md_simulator::Basic::add_force_field)
    .def("add_constraint",&md_simulator::Basic::add_constraint)
    .def("add_writer",&md_simulator::Basic::add_writer)

    .def_readwrite("time",&md_simulator::Basic::time)      
    .def_readwrite("dt",&md_simulator::Basic::dt) 
    .def_readwrite("initial_step",&md_simulator::Basic::initial_step) 
    .def_readwrite("final_step",&md_simulator::Basic::final_step) 
    .def_readwrite("initial_time",&md_simulator::Basic::initial_time) 
    .def_readwrite("final_time",&md_simulator::Basic::final_time) 
    .def_readwrite("initialized",&md_simulator::Basic::initialized) 
    .def_readwrite("use_time",&md_simulator::Basic::use_time)      
    .def_readwrite("use_step",&md_simulator::Basic::use_step)      

    .add_property("atom_data", &md_simulator::Basic::get_atom_data, &md_simulator::Basic::set_atom_data)
    .add_property("integrator", &md_simulator::Basic::get_integrator, &md_simulator::Basic::set_integrator)
  ;

}



} //md_simulator
} //objects
} // namespace caviar

