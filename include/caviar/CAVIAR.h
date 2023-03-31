
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

#ifndef CAVIAR_CAVIAR_H
#define CAVIAR_CAVIAR_H

#include "caviar/utility/caviar_config.h"
#include "caviar/utility/types.h"

#include <iostream>
#include <fstream>
#include <map>

#if defined(CAVIAR_WITH_MPI)
#include <mpi.h>
#endif

#include "caviar/utility/common_template_functions.h"

CAVIAR_NAMESPACE_OPEN
namespace interpreter {
class Communicator;
class Error; 
class Output; 
class Input;
class Object_handler;
class Object_container;
class Object_creator;
}

/**
 * This class is the base of CAVIAR package. It contains all the input/output
 * classes and also the related class of objects. One instance of this class is 
 * created at the begining of any CAVIAR calculations.
 */
class CAVIAR {
public:

#if defined(CAVIAR_WITH_MPI)
  /**
   * Constructor in mpi mode
   */
  CAVIAR (int, char**, MPI_Comm);
#else

  /**
   * Constructor in serial mode
   */
  CAVIAR (int, char**);
#endif

  /**
   * Destructor.
   */
  ~CAVIAR ();

  /**
   *  the function to call all the calculations.
   */
  void execute ();
  
#if defined(CAVIAR_WITH_MPI)
  MPI_Comm mpi_comm;
#endif
  class interpreter::Communicator *comm;
  class interpreter::Error *error;
  class interpreter::Output *output;
  class interpreter::Input *input;
  class interpreter::Object_handler *object_handler;
  class interpreter::Object_container *object_container;
  class interpreter::Object_creator *object_creator;
  
  std::ofstream log;
  std::istream in;
  std::ostream out, err;
  bool log_flag, out_flag, err_flag;
  int argc;
  char **argv;

  // these are some helper variables for the interpreter. Not for users.
  // They can be public.
  int interpreter_num_Input_class;
  bool interpreter_break_called;
  bool interpreter_continue_called;
};

CAVIAR_NAMESPACE_CLOSE

#endif
