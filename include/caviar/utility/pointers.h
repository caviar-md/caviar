
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

#ifndef CAVIAR_POINTERS_H
#define CAVIAR_POINTERS_H

#ifdef CAVIAR_WITH_MPI
#include <mpi.h>
#endif

#include "caviar/CAVIAR.h"

namespace caviar {
namespace interpreter {
class Parser;
class Lexer;
class Communicator;
class Error; 
class Output; 
class Input;
class Object_handler;
class Object_container;
class Object_creator;
}

// note that this class only have pointer type objects. If you want to add anything
// to all the base objects, add them to the FC_BASE_OBJECT_COMMON_TOOLS definition
// in the 'macro_function.h'
class Pointers {
public:
  inline constexpr Pointers (class CAVIAR *fptr) :
    fptr {fptr},
#if defined(CAVIAR_WITH_MPI)
    mpi_comm {fptr->mpi_comm},
#endif
    comm {fptr->comm},
    error {fptr->error},
    output {fptr->output},
    input {fptr->input},
    object_handler {fptr->object_handler},
    object_container {fptr->object_container},
    object_creator {fptr->object_creator},    
    log {fptr->log},
    in {fptr->in},
    out {fptr->out},
    err {fptr->err},
    log_flag {fptr->log_flag},
    out_flag {fptr->out_flag},
    err_flag {fptr->err_flag} {}

protected:
  inline ~Pointers () = default;  // destructor need not be virtual b/c it is public
public:
  class CAVIAR *fptr;
#if defined(CAVIAR_WITH_MPI)
  MPI_Comm &mpi_comm;
#endif
  class interpreter::Communicator *&comm;
  class interpreter::Error *&error;
  class interpreter::Output *&output;
  class interpreter::Input *&input;
  class interpreter::Object_handler *&object_handler;
  class interpreter::Object_container *&object_container;
  class interpreter::Object_creator *&object_creator;
  
  std::ofstream &log;
  std::istream &in;
  std::ostream &out, &err;
  bool &log_flag, &out_flag, &err_flag;
};

} // namespace caviar

#endif
