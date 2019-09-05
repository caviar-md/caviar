
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

#ifndef CAVIAR_INTERPRETER_INPUT_H
#define CAVIAR_INTERPRETER_INPUT_H

#include "caviar/utility/pointers.h"
#include <istream>
#include <map>

namespace caviar {

// a pointer to boolean function of Input class.
//using InputCommandFunc = bool (Input::*) (class caviar::interpreter::Parser *); 
namespace interpreter {
using InputCommandFunc = char (Input::*) (class caviar::interpreter::Parser *); 

/**
 * This class contains input streams lexer and parser and handles the 
 * basic commands of CAVIAR scripting languages.
 * 
 */
class Input : public Pointers {
public:

  /**
   * Constructor.
   */
  Input (class CAVIAR *);

  /**
   * Constructor with a input file name @p filename
   */
  Input (class CAVIAR *, const std::string & filename);

  /**
   * Constructor with a string stream
   */
  Input (class CAVIAR *, std::istringstream & stream);
  ~Input ();
  
  /**
   *  the function to start reading
   */
  void read ();
public:
  class caviar::interpreter::Parser *parser;
  class CAVIAR *fptr;

  const static std::map<std::string,InputCommandFunc> commands_map;

  bool read (caviar::interpreter::Parser *);

  bool read_command (Parser *);  

  // commands


  char command_output (Parser *);
  char command_object_container (Parser *);
  char command_read_script_from_file (Parser *);
  char command_exit_program (Parser *);
  char command_echo(Parser *);
  char command_print(Parser *);
  char command_fprint(Parser *);

  char command_for(Parser *);
  char command_next(Parser *);

  char command_while(Parser *);
  char command_do(Parser *);
  char command_end_do(Parser *);

  char command_break(Parser *);
  char command_continue(Parser *);

  char command_if(Parser *);
  char command_else_if(Parser *);
  char command_else(Parser *);
  char command_end_if(Parser *);

  char command_evaluate(Parser *);
  char command_evaluate(const std::string &);

  char command_calculate(Parser *);
  char command_calculate(const std::string &);

  char command_compare_string (Parser *);
  char command_compare_string (const std::string &);
  char command_compare_real (Parser *);
  char command_compare_real (const std::string &);
  char command_compare_int (Parser *);
  char command_compare_int (const std::string &);
  char command_compare (Parser *);
  char command_compare (const std::string &);

  char command_delete_object(Parser *);
  char command_delete_object(const std::string &);

  char command_help(Parser *);
  char command_exit(Parser *);
  char command_delete(Parser *);

  char command_include(Parser *);
  char command_import(Parser *);
  char command_read (caviar::interpreter::Parser *);


  char command_function(Parser *);
  char command_end_function(Parser *);

  char command_class(Parser *);
  char command_end_class(Parser *);

};
} //interpreter
} // namespace caviar

#endif

