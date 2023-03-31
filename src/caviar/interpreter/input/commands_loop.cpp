
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
#include "caviar/interpreter/object_handler.h"
#include "caviar/utility/interpreter_io_headers.h"

#include <map>
#include <cmath>
#include <algorithm>

CAVIAR_NAMESPACE_OPEN
namespace interpreter
{
  // used in order to remove 'if' or 'elseif' from the line and get the condition
  static std::string remove_the_first_word(const std::string &str)
  {
    int i = 0;
    while (isblank(str[i]))
      i++;
    while (!isblank(str[i]))
      i++;
    while (isblank(str[i]))
      i++;
    if (i > static_cast<int>(str.size()))
      return "";
    return str.substr(i);
  }
  // commands

  char Input::command_do(Parser *parser)
  {
    //----------------
    // loop interpreter
    //----------------
    std::string loop_initial_condition;
    std::string loop_final_condition;
    std::string loop_text;

    int no_nested_loop = 0;

    int condition_type = 0; // 0: no_condition; -1: initial_condition; +1:final_condition
    //-----------
    // read loop
    //-----------

    loop_initial_condition = remove_the_first_word(parser->line);
    if (!parser->go_to_next_line())
      error->all(FC_FILE_LINE_FUNC_PARSE, "expected 'enddo' but reached 'eof'.");

    while (true)
    {
      auto t = parser->get_raw_token();
      auto ts = t.string_value;
      if (t.kind == caviar::interpreter::Kind::eof)
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "expected 'enddo' but reached 'eof'.");
      }
      else if (t.kind == caviar::interpreter::Kind::eol)
      {
        continue;
      }
      else if (t.kind == caviar::interpreter::Kind::identifier)
      {
        if (ts == "do")
        {
          ++no_nested_loop;
        }
        else if (ts == "enddo" || ts == "end_do")
        {
          if (no_nested_loop > 0)
          {
            --no_nested_loop;
          }
          else
          {
            loop_final_condition = remove_the_first_word(parser->line);
            break;
          }
        }
      }
      loop_text += parser->line;
      if (!parser->go_to_next_line())
        error->all(FC_FILE_LINE_FUNC_PARSE, "expected 'enddo' but reached 'eof'.");
    }

    //---------
    // print loop
    //---------

    // std::cout << "loop_initial_condition: " << loop_initial_condition << std::endl;
    // std::cout << "loop_text:\n" << loop_text << std::endl;
    // std::cout << "loop_final_condition: " << loop_final_condition << std::endl;

    //--------
    // run loop
    //--------

    if (loop_initial_condition != "")
      condition_type = -1;
    // if (loop_final_condition!="") condition_type = +1;
    // std::cout << "condition_type: " << condition_type << std::endl;
    // if (loop_initial_condition!="" && loop_final_condition != "")
    // error->all(FC_FILE_LINE_FUNC_PARSE, "having two conditions in one loop is not supported.");

    while (true)
    {
      if (condition_type == -1)
        if (!command_compare_real(loop_initial_condition))
          break;

      command_evaluate(loop_text);
      if (fptr->interpreter_break_called)
      {
        fptr->interpreter_break_called = false;
        break;
      }
      if (fptr->interpreter_continue_called)
      {
        fptr->interpreter_continue_called = false;
        continue;
      }
      // if (condition_type == +1)
      // if (!compare_real(loop_final_condition)) break;
    }

    return true;
  }

  char Input::command_end_do(Parser *parser)
  {
    error->all(FC_FILE_LINE_FUNC_PARSE, "having 'end_do' without 'do' is not permitted ");
    return true;
  }

  char Input::command_while(Parser *)
  {
    error->all(FC_FILE_LINE_FUNC_PARSE, "having 'while' without 'do' is not permitted ");
    return true;
  }

  char Input::command_break(Parser *)
  {
    fptr->interpreter_break_called = true;
    return false;
  }

  char Input::command_continue(Parser *)
  {
    fptr->interpreter_continue_called = true;
    return false;
  }

} // interpreter
CAVIAR_NAMESPACE_CLOSE
