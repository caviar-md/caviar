
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

//#define DEBUG_ME

CAVIAR_NAMESPACE_OPEN
namespace interpreter {
// used in order to remove 'if' or 'elseif' from the line and get the condition
static std::string remove_the_first_word (const std::string &str) {
  int i = 0;
  while (isblank(str[i])) i++;
  while (!isblank(str[i])) i++;
  while (isblank(str[i])) i++;
  if (i > static_cast<int>(str.size())) return "";
  return str.substr(i);
}
// commands


char Input::command_if (Parser *parser) {
  //-------------
  // if_interpreter
  //-------------
  bool if_part = true;
  std::string if_condition;
  std::string if_text;

  bool elseif_part = false;
  int elseif_index = -1;
  // there can be more than one 'elseif'
  std::vector<std::string> elseif_condition;
  std::vector<std::string> elseif_text;

  bool else_found = false;
  bool else_part = false;;
  std::string else_text;

  int no_nested_if = 0;

  //--------
  // read if
  //--------

  if_condition = remove_the_first_word(parser->line) ;
  if (!parser->go_to_next_line())
    error->all(FC_FILE_LINE_FUNC_PARSE, "expected 'endif' but reached 'eof'.");

  while (true) {

    auto t = parser -> get_raw_token();
    auto ts = t.string_value;

    if (t.kind == caviar::interpreter::Kind::eof) {
      error->all(FC_FILE_LINE_FUNC_PARSE, "expected 'endif' but reached 'eof'.");
    } else if (t.kind == caviar::interpreter::Kind::eol) {
      continue;
    } else if (t.kind == caviar::interpreter::Kind::identifier) {
      if (ts == "if") {
        ++no_nested_if; 
      } else if (ts == "endif" || ts == "end_if") {
        if (no_nested_if > 0) {
          --no_nested_if;
        } else {
          break;
        }
      } else if (ts == "else") {
        if (no_nested_if > 0) {

        } else {
          if_part = false;
          elseif_part = false;
          else_found = true;
          else_part = true;
          if (!parser->go_to_next_line())
            error->all(FC_FILE_LINE_FUNC_PARSE, "expected 'endif' but reached 'eof'.");
          continue;

        }
      } else if (ts == "elseif" || ts == "else_if") {
        if (no_nested_if > 0) {

        } else {

          if (else_found)
            error->all(FC_FILE_LINE_FUNC_PARSE, "can not have 'elseif' after 'else'.");
          if_part = false;
          elseif_part = true;
          elseif_condition.emplace_back(remove_the_first_word(parser->line)) ;
          elseif_text.emplace_back("");
          elseif_index++;
          if (!parser->go_to_next_line())
            error->all(FC_FILE_LINE_FUNC_PARSE, "expected 'endif' but reached 'eof'.");
          continue;

        }
      }
    }

    if (if_part) {
      if_text += parser->line ;
    } else if (elseif_part) {
      elseif_text[elseif_index] += parser->line;
    } else if (else_part) {
      else_text += parser->line;
    }

    if (!parser->go_to_next_line())
      error->all(FC_FILE_LINE_FUNC_PARSE, "expected 'endif' but reached 'eof'.");
  }

  //---------
  // print_if
  //---------
  /*
  std::cout << "if_condition: " << if_condition << std::endl;
  std::cout << "if_text:\n" << if_text << std::endl;
  for (int i = 0; i < elseif_condition.size(); ++i) {
    std::cout << "elseif_condition: " << elseif_condition[i] << std::endl;
    std::cout << "elseif_text:\n" << elseif_text[i] << std::endl;
  }
  std::cout << "else_text:\n" << else_text << std::endl;
  */
  //--------
  // run_if
  //--------

  if (command_compare_real(if_condition)) {
    command_evaluate(if_text);
    return true;
  }

  for (unsigned int i = 0; i < elseif_condition.size(); ++i) {
    if (command_compare_real(elseif_condition[i])) {
      command_evaluate(elseif_text[i]);
      return true;
    }
  }

  if (else_found) {
    command_evaluate(else_text);
    return true;
  }

  return true;
}



char Input::command_else_if (Parser *) {
  error->all (FC_FILE_LINE_FUNC_PARSE,"catched 'elseif' without 'if'");
  return true;
}

char Input::command_else (Parser *) {
  error->all (FC_FILE_LINE_FUNC_PARSE,"catched 'else' without 'if'");
  return true;
}

char Input::command_end_if (Parser *) {
  error->all (FC_FILE_LINE_FUNC_PARSE,"catched 'endif' without 'if'");
  return true;
}


char Input::command_for(Parser*) {
  return true;
}
char Input::command_next(Parser*) {
  return true;
}


} //interpreter
CAVIAR_NAMESPACE_CLOSE

