
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

namespace caviar {
namespace interpreter {
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
// commands
*/
char Input::command_read_script_from_file (Parser *parser) {
  auto file_name = parser->get_string();
  std::string st = "read_script_from_file : " +  file_name;
  output->info(st);
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: Input::command_read_script_from_file : create INPUT" << std::endl;
#endif
  caviar::interpreter::Input inp (fptr, file_name);
  inp.read();
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: Input::command_read_script_from_file : delete INPUT" << std::endl;
#endif
  return true;
}

char Input::command_help (Parser *parser) {
  parser->end_of_line();
  std::cout << "Science Helps People!\n";
  return true;
}

char Input::command_exit_program (Parser *) {
  return false;
}

char Input::command_output (Parser *parser) {
  return output->read (parser);
}

char Input::command_object_container (Parser *) {
  return object_container->read (parser);
}

char Input::command_echo (Parser * parser) {
  std::string st = "";

  while (true) {
    //auto t = parser -> get_val_token();
    auto t = parser -> get_raw_token();
    if (t.kind == caviar::interpreter::Kind::eof) break;
    else if (t.kind == caviar::interpreter::Kind::eol) break;
    else if (t.kind == caviar::interpreter::Kind::string) {
      // having 'get_string()' here may complicate 'kind::plus' conditions
      st += t.string_value;
    } else if (t.kind == caviar::interpreter::Kind::identifier) {
      auto var_name = t.string_value;
      std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;
      it = object_container->dictionary.find(var_name);
      if (it == object_container->dictionary.end())   {
        error->all (FC_FILE_LINE_FUNC_PARSE, "unknown object name :" + var_name);      
      } else if (it->second.type == object_handler::gdst( "boolean_variable" )) {
        parser-> keep_current_token();
        bool r = parser-> get_bool();
        st += std::to_string(r);
      } else if (it->second.type == object_handler::gdst( "int_variable" )) {
        parser-> keep_current_token();
        double r = parser-> get_real();
        st += std::to_string(r);
      } else if (it->second.type == object_handler::gdst( "real_variable" )) {
        parser-> keep_current_token();
        double r = parser-> get_real();
        st += std::to_string(r);
      } else if (it->second.type == object_handler::gdst( "real_3d_vector" )) {
        auto v = object_container -> real_3d_vector [it->second.index];        
        st += "{" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", "
            +  std::to_string(v.z) + "}";
      } else if (it->second.type == object_handler::gdst( "int_3d_vector" )) {
        auto v = object_container -> int_3d_vector [it->second.index];        
        st += "{" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", "
            +  std::to_string(v.z) + "}";
      } else if (it->second.type == object_handler::gdst( "string_variable" )) {
        //st += object_container->string_variable[it->second.index];
        parser-> keep_current_token();
        std::string r = parser-> get_string();
        st += r;
      } else {
        error->all (FC_FILE_LINE_FUNC_PARSE, "this echo type is not implemented yet");      
      }
    } else if (t.kind == caviar::interpreter::Kind::real_number || t.kind == caviar::interpreter::Kind::int_number) {
      // the distinction between 'int_number' and 'real_number' is not nessecary
      // and it just complicates the code.
      parser-> keep_current_token();
      double r = parser-> get_real();
      st += std::to_string(r);
    } else if (t.kind == caviar::interpreter::Kind::plus || t.kind == caviar::interpreter::Kind::minus) {
      // the same as above
      parser-> keep_current_token();
      double r = parser-> get_real();
      st += std::to_string(r);
    } else {
      error->all (FC_FILE_LINE_FUNC_PARSE, "this echo type is not implemented yet");      
    }
  }

  std::cout << st << std::endl;
  return true;
}

char Input::command_print (Parser *parser) {
  std::string st = parser -> get_string();
  std::cout << st << std::endl;

//  std::string dollars = "$";
//  std::size_t found = st.find(dollars);

  return true;
}

char Input::command_evaluate(const std::string &str) {
  std::string st = str + "\n";
  std::istringstream iss (st);
  //std::cout << "iss: " << iss.str() << std::endl;
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: Input::command_evaluate : create INPUT" << std::endl;
#endif
  caviar::interpreter::Input inp (fptr, iss);
  inp.read();
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: Input::command_evaluate : delete INPUT" << std::endl;
#endif
  return true;
}

char Input::command_evaluate (Parser *parser) {
  command_evaluate ( parser -> get_string() );
  return true;
}

// TODO: complete
char Input::command_compare_int (Parser *parser) {
  return command_compare_real(parser);
}

char Input::command_compare_real (const std::string &st) {
  std::istringstream iss (st);
  Parser p(fptr, iss);
  //std::cout << "X: " << st << " , result: " << p.compare_real() << std::endl;
  return p.compare_real();
}

char Input::command_compare_real (Parser *parser) {
  std::string st = parser->rest_of_line();
  std::istringstream iss (st);
  Parser p(fptr, iss);
  bool result =  p.compare_real();
  //std::cout << "X: " << st << " , result: " << result << std::endl;
  parser->go_to_next_line();
  return result;
}

// TODO: complete
char Input::command_compare (Parser *parser) {
  return command_compare_real(parser);
}

char Input::command_calculate(const std::string &str) {
  std::string st = str + "\n";
  std::istringstream iss (st);
  //std::cout << "iss: " << iss.str() << std::endl;
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: Input::command_calculate : create INPUT" << std::endl;
#endif
  caviar::interpreter::Input inp (fptr, iss);
  inp.read();
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: Input::command_calculate : delete INPUT" << std::endl;
#endif
  return true;
}

char Input::command_calculate (Parser *parser) {
  int num_of_condition_operators = 0;
  int condition_operator = 10;
  double left_side = 0;
  double right_side = 0;
  bool operator_called = false;
  while (true) {
    auto t = parser -> get_raw_token();
    if (t.kind == caviar::interpreter::Kind::eof) break;
    else if (t.kind == caviar::interpreter::Kind::eol) break;
    else if (t.kind == caviar::interpreter::Kind::smaller) { //-2
      ++num_of_condition_operators; operator_called = true;
      condition_operator = -2;
    } else if (t.kind == caviar::interpreter::Kind::larger) { //+2
      ++num_of_condition_operators; operator_called = true;
      condition_operator = +2;
    } else if (t.kind == caviar::interpreter::Kind::eqsmaller) { //-1
      ++num_of_condition_operators; operator_called = true;
      condition_operator = -1;
    } else if (t.kind == caviar::interpreter::Kind::eqlarger) { //+1
      ++num_of_condition_operators; operator_called = true;
      condition_operator = +1;
    } else if (t.kind == caviar::interpreter::Kind::equal) { // 0
      ++num_of_condition_operators; operator_called = true;
      condition_operator = 0;
    } else {
      parser -> keep_current_token();
      if (operator_called) {right_side = parser->get_real(); break;}
      else left_side = parser->get_real();
    } 
  }

  if (num_of_condition_operators == 0)
    error->all (FC_FILE_LINE_FUNC_PARSE, "expected a conditional operator" ); 

  if (num_of_condition_operators > 1)
    error->all (FC_FILE_LINE_FUNC_PARSE, "only one conditional operator is allowed" ); 

  //std::cout <<"left: " << left_side << " , right: " << right_side << std::endl;

  if (condition_operator == -2) return (left_side <  right_side);
  if (condition_operator == -1) return (left_side <= right_side);
  if (condition_operator ==  0) return (left_side == right_side);
  if (condition_operator ==  1) return (left_side >= right_side);
  if (condition_operator ==  2) return (left_side >  right_side);

  return true;
}


char Input::command_function(Parser*) {
  return true;
}
char Input::command_end_function(Parser*) {
  return true;
}
char Input::command_class(Parser*) {
  return true;
}
char Input::command_end_class(Parser*) {
  return true;
}
char Input::command_read (caviar::interpreter::Parser*) {
  return true;
}
char Input::command_include(Parser*) {
  return true;
}
char Input::command_import(Parser*) {
  return true;
}
char Input::command_exit(Parser*) {
  //exit(0);
  return false;
}


char Input::command_compare_string(Parser*) {
  return true;
}
char Input::command_delete(Parser*) {
  return true;
}

} //interpreter
} // namespace caviar

