
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

#include "caviar/CAVIAR.h"

#include <string>

#include "caviar/interpreter/all.h"

namespace caviar {

#if defined(CAVIAR_WITH_MPI)
CAVIAR::CAVIAR (int argc, char **argv, MPI_Comm mpi_comm) :
  mpi_comm {mpi_comm},
#else
CAVIAR::CAVIAR (int argc, char **argv) :
#endif

  comm {new interpreter::Communicator {this}},
  error {new interpreter::Error {this}},
  output {new interpreter::Output {this}},
  input {new interpreter::Input {this}},
  object_handler {new interpreter::Object_handler {this}},
  object_container {new interpreter::Object_container {this}},
  object_creator {new interpreter::Object_creator {this}},  
  in {std::cin.rdbuf()},
  out {std::cout.rdbuf()},
  err {std::cerr.rdbuf()},
  log_flag {true},
  out_flag {true},
  err_flag {true},
  argc{argc},
  argv{argv} {
    if (comm->me == 0) log.open ("log");

    // this value is set to '1' because one input class is already constructed
    interpreter_num_Input_class = 1;

    interpreter_break_called = false;
    interpreter_continue_called = false;
}

#if defined(CAVIAR_WITH_MPI)

#else
CAVIAR::CAVIAR (std::string str) :

  comm {new interpreter::Communicator {this}},
  error {new interpreter::Error {this}},
  output {new interpreter::Output {this}},
  input {new interpreter::Input {this}},
  object_handler {new interpreter::Object_handler {this}},
  object_container {new interpreter::Object_container {this}},
  object_creator {new interpreter::Object_creator {this}},  
  in {std::cin.rdbuf()},
  out {std::cout.rdbuf()},
  err {std::cerr.rdbuf()},
  log_flag {true},
  out_flag {true},
  err_flag {true},
  construct_str{str}
  {
    if (comm->me == 0) log.open ("log");

    // this value is set to '1' because one input class is already constructed
    interpreter_num_Input_class = 1;

    interpreter_break_called = false;
    interpreter_continue_called = false;
}
#endif

CAVIAR::~CAVIAR () {
  delete input;
  delete output;
  delete error;
  delete comm;
  delete object_handler;  
  delete object_container;    
  delete object_creator;    
}

void CAVIAR::execute () {
  std::string greeting = "CAVIAR-";
  greeting += std::to_string (CAVIAR_MAJOR_VERSION);
  greeting += ".";
  greeting += std::to_string (CAVIAR_MINOR_VERSION);
  greeting += ".";
  greeting += std::to_string (CAVIAR_PATCH_VERSION);
    
  output->info(greeting);
    
  input->read ();
}


} // namespace caviar

#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(_caviar)
{
    class_<caviar::CAVIAR,boost::noncopyable>("caviar", init<std::string>())
        .def("greet", &caviar::CAVIAR::execute)

    ;
};

