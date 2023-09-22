
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

#include "caviar/interpreter/input.h"
#include "caviar/interpreter/object_creator.h"
#include "caviar/interpreter/input/commands_map.h"
#include "caviar/interpreter/object_handler.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/utility/file_utility.h"

#include <map>
#include <cmath>
#include <algorithm>
#include <cstring>

//#define DEBUG_ME

CAVIAR_NAMESPACE_OPEN
namespace interpreter
{
  // used in order to remove 'if' or 'elseif' from the line and get the condition
  /*
  static std::string remove_the_first_word (const std::string &str) {
    int i = 0;
    while (isblank(str[i])) i++;
    while (!isblank(str[i])) i++;
    while (isblank(str[i])) i++;
    if (i > static_cast<int>(str.size())) return "";
    return str.substr(i);
  }
  */
  //0     1     2
  //cav  -i  /dsaads/das
  Input::Input(CAVIAR *fptr, int argc, char **argv) : Pointers{fptr}, fptr{fptr}
  {
    fptr->interpreter_num_Input_class++;

    bool input_file_found = false;
    std::string input_file;

    for (int i = 0; i < argc; ++i)
    {
      //std::cout << "argv " << i << " " << argv[i] << " " << std::endl;
      if (std::strcmp(argv[i],"-i") == 0)
      {
        if (argc > i + 1)
        {
          input_file = argv[i+1];
          if (file_exists_1(input_file))
          {
            input_file_found = true;
            input_file_directory = directory_of_file(input_file);
          }
          else
          {
            std::cout << "Error: File doesn't exist " << input_file << std::endl;
            exit(1);
          }
        }
        else
        {
          std::cout << "Error: Expected input file after '-i' argument" << std::endl;
          exit(1);
        }
      }
    }

    if (input_file_found)
      parser = new Parser{fptr, input_file};
    else
      parser = new Parser{fptr};

#ifdef DEBUG_ME
    std::cout << "DEBUG_ME: INPUT constructor t0" << std::endl;
#endif
  }

  Input::Input(CAVIAR *fptr) : Pointers{fptr}, parser{new Parser{fptr}},
                               fptr{fptr}
  {
    fptr->interpreter_num_Input_class++;
#ifdef DEBUG_ME
    std::cout << "DEBUG_ME: INPUT constructor t1" << std::endl;
#endif
  }

  Input::Input(CAVIAR *fptr, const std::string &file) : Pointers{fptr},
                                                        parser{new Parser{fptr, file}}, fptr{fptr}
  {
    fptr->interpreter_num_Input_class++;
#ifdef DEBUG_ME
    std::cout << "DEBUG_ME: INPUT constructor t2" << std::endl;
#endif
  }

  Input::Input(CAVIAR *fptr, std::istringstream &iss) : Pointers{fptr},
                                                        parser{new Parser{fptr, iss}}, fptr{fptr}
  {
    fptr->interpreter_num_Input_class++;
#ifdef DEBUG_ME
    std::cout << "DEBUG_ME: INPUT constructor t3" << std::endl;
#endif
  }

  Input::~Input()
  {
    fptr->interpreter_num_Input_class--;
#ifdef DEBUG_ME
    std::cout << "DEBUG_ME: INPUT destructor " << std::endl;
#endif
    delete parser;
  }

  // called by CAVIAR object at execute()
  void Input::read()
  {
    while (read_command(parser))
      ;
    if (fptr->interpreter_num_Input_class == 1)
    {
      if (fptr->interpreter_break_called)
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "stray 'break' command.");
      }
      else if (fptr->interpreter_continue_called)
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "stray 'continue' command.");
      }
    }
  }

  bool Input::read(caviar::interpreter::Parser *parser)
  {
    while (read_command(parser))
      ;
    return true;
  }

  bool Input::read_command(Parser *parser)
  {
    // there can be more than one Input levels when there's a 'break' or 'continue'
    // call. For example in cases there's an 'if' condition inside a 'do' loop.
    // In that case, the interpreter isn't allowed to do anything unless the command
    // is handled.
    if (fptr->interpreter_break_called || fptr->interpreter_continue_called)
    {
      return false;
    }
    Token t;
    do
    { // go to the first non-eol token
      t = parser->get_raw_token();
    } while (t.kind == caviar::interpreter::Kind::eol);

    std::string command;
    if (t.kind == caviar::interpreter::Kind::eof)
      return false;
    if (t.kind == caviar::interpreter::Kind::identifier)
      command = t.string_value;
    auto command_lowercase = command;
#ifdef CAVIAR_SCRIPT_COMMAND_CASE_INSENSITIVE
    // transform 'command' to lower case
    std::transform(command_lowercase.begin(), command_lowercase.end(),
                   command_lowercase.begin(), ::tolower);
#endif
    if (commands_map.count(command_lowercase) != 0)
      return (this->*commands_map.at(command_lowercase))(parser);
    else if (object_creator->commands_map.count(command_lowercase) != 0)
    {
      return (object_creator->*object_creator->commands_map.at(command_lowercase))(parser);
    }
    else if (object_handler->commands_map.count(command_lowercase) != 0)
    {
      return (object_handler->*object_handler->commands_map.at(command_lowercase))(parser);
    }
    else if (object_container->all_names.count(command) != 0)
    { // object name is case sensitive
      return object_handler->read_object(parser, command);
    }
    else
    {
      error->all(FC_FILE_LINE_FUNC_PARSE, static_cast<std::string>("Invalid command or object name: ") + command_lowercase);
    }
    return true;
  }
} // interpreter

CAVIAR_NAMESPACE_CLOSE
