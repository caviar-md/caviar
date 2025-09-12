
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

#ifndef CAVIAR_INTERPRETER_INPUT_COMMANDSMAP_H
#define CAVIAR_INTERPRETER_INPUT_COMMANDSMAP_H

#include "caviar/interpreter/input.h"

CAVIAR_NAMESPACE_OPEN
namespace interpreter {
//using CommandFunc = bool (Input::*) (class caviar::interpreter::Parser *); 

/**
 * A map between command names and the related functions.
 * 
 * 
 */
const std::map<std::string,InputCommandFunc> Input::commands_map = {

  {"object_container", &Input::command_object_container},

  {"function",      &Input::command_function},
  {"end_function",  &Input::command_end_function},

  {"class",        &Input::command_class},
  {"end_class",    &Input::command_end_class},

  {"read",          &Input::command_read},

  {"include",       &Input::command_include},

  {"import",        &Input::command_import},     

  {"output",        &Input::command_output},

  {"exit",          &Input::command_exit},
  {"quit",          &Input::command_exit},

  {"print",         &Input::command_print},
  {"fprint",         &Input::command_print},


  {"echo",          &Input::command_echo},

  {"for",           &Input::command_for},
  {"next",          &Input::command_next},

  {"do",            &Input::command_do},
  {"end_do",        &Input::command_end_do},  
  {"enddo",        &Input::command_end_do},  

  {"while",         &Input::command_while},

  {"break",         &Input::command_break},
  {"continue",      &Input::command_continue},

  {"if",            &Input::command_if},
  {"endif",         &Input::command_end_if},
  {"end_if",        &Input::command_end_if},
  {"else_if",       &Input::command_else_if},
  {"elseif",        &Input::command_else_if},
  {"else",          &Input::command_else},

  {"evaluate",      &Input::command_evaluate},
  {"calculate",     &Input::command_calculate},

  {"compare_real",  &Input::command_compare_real},
  {"compare_int",   &Input::command_compare_int},
  {"compare_string",&Input::command_compare_string},
  {"compare",       &Input::command_compare},

  {"delete",        &Input::command_delete},

  {"help",          &Input::command_help},
};

} //interpreter
CAVIAR_NAMESPACE_CLOSE

#endif
