
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

#include "caviar/interpreter/input/parser.h"
#include "caviar/interpreter/input/lexer.h"
#include "caviar/interpreter/error.h"
#include "caviar/utility/interpreter_io_headers.h"

#include <cmath>

//#define DEBUG_ME

namespace caviar {
namespace interpreter {
// this operator is used in debugging parser.
std::ostream& operator<<(std::ostream& out, const Kind value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
      PROCESS_VAL(Kind::identifier); PROCESS_VAL(Kind::string); PROCESS_VAL(Kind::int_number);
      PROCESS_VAL(Kind::real_number); PROCESS_VAL(Kind::eof); PROCESS_VAL(Kind::plus);
      PROCESS_VAL(Kind::minus); PROCESS_VAL(Kind::mul); PROCESS_VAL(Kind::div);
      PROCESS_VAL(Kind::truediv); PROCESS_VAL(Kind::modulus); 
      PROCESS_VAL(Kind::assign); PROCESS_VAL(Kind::asplus); PROCESS_VAL(Kind::asminus);
      PROCESS_VAL(Kind::asmul); PROCESS_VAL(Kind::asdiv);
      PROCESS_VAL(Kind::lp); PROCESS_VAL(Kind::rp); PROCESS_VAL(Kind::squote); PROCESS_VAL(Kind::dquote);
      PROCESS_VAL(Kind::comment); PROCESS_VAL(Kind::backslash); PROCESS_VAL(Kind::eol);
      PROCESS_VAL(Kind::lcb); PROCESS_VAL(Kind::rcb);
      PROCESS_VAL(Kind::lsb); PROCESS_VAL(Kind::rsb);
      PROCESS_VAL(Kind::eqsmaller); PROCESS_VAL(Kind::eqlarger);  
      PROCESS_VAL(Kind::equal);    PROCESS_VAL(Kind::smaller); PROCESS_VAL(Kind::larger);
      PROCESS_VAL(Kind::unknown);
      default: s = "Kind::unknown"; break;
    }
#undef PROCESS_VAL

    return out << s;
}

Parser::Parser (CAVIAR *fptr) : Pointers{fptr}, token_stream{new Token_stream{fptr}}, get_new_token{true}, line{token_stream->line}, col{token_stream->col} {
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: PARSER constructor t1" << std::endl;
#endif
}

Parser::Parser (CAVIAR *fptr, const std::string &file) : Pointers{fptr}, token_stream{new Token_stream{fptr, file}}, get_new_token{true}, line{token_stream->line}, col{token_stream->col} {
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: PARSER constructor t2" << std::endl;
#endif
}

Parser::Parser (CAVIAR *fptr, std::istringstream &iss) : Pointers{fptr}, token_stream{new Token_stream{fptr, iss}}, get_new_token{true}, line{token_stream->line}, col{token_stream->col} {
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: PARSER constructor t3" << std::endl;
#endif
}

Parser::~Parser () {
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: PARSER destructor" << std::endl;
#endif
  delete token_stream;
}

bool Parser::end_of_line () {
  Token token = get_raw_token();
  return token.kind == caviar::interpreter::Kind::eol;
}

bool Parser::assignment () {
  Token token = get_raw_token();
  return token.kind == caviar::interpreter::Kind::assign;
}

Token Parser::get_raw_token () {
  if (get_new_token) {
/*
    while (true) {
      auto t = token_stream->get();
      if (t.kind == caviar::interpreter::Kind::comment) {
        //std::cout << "X: comment found at line: " << line << std::endl;
        continue;
      }
      return t;
    }
*/
   return token_stream->get();
  
 } else {
    get_new_token = true;
    return token_stream->current();
  }
}

Token Parser::get_val_token () {
  auto token = get_raw_token ();
  switch (token.kind) {
    case Kind::identifier:
    case Kind::string:
    case Kind::int_number:
    case Kind::real_number:
    case Kind::eol:
    case Kind::eof:
      return token;
    
    case Kind::plus:
      token = get_raw_token ();
      switch (token.kind) {
        case Kind::int_number:
        case Kind::real_number:
          return token;
        default:
          error->all (FC_FILE_LINE_FUNC_LINE_COL, "Unexpected token");
      }
    
    case Kind::minus:
      token = get_raw_token ();
      switch (token.kind) {
        case Kind::int_number:
          token.int_value *= -1;
          break;
        case Kind::real_number:
          token.real_value *= -1;
          break;
        default:
          error->all (FC_FILE_LINE_FUNC_LINE_COL, "Unexpected token");
      }
      return token;
    

    default:
      std::cout << token.kind << std::endl; 
      error->all (FC_FILE_LINE_FUNC_LINE_COL, "Unexpected token");
      return token; // WARNING      
  }

}


std::ostream & Parser::write_to_stream (std::ostream &stream) {
  auto token = get_raw_token();
  if (token.kind == caviar::interpreter::Kind::string) stream << token.string_value << std::endl;
  else if (token.kind == caviar::interpreter::Kind::identifier) {
    //if (string_variables.count (token.string_value))
    //stream << string_variables.at(token.string_value) << std::endl;//XXX
  } else stream << expression(false) << std::endl;
  return stream;
}

std::string Parser::get_command_identifier () {
/*
  Token token;
  do { // go to the first non-eol token
    token = get_raw_token();
    if (token.kind == caviar::interpreter::Kind::comment) continue;
  } while (token.kind == caviar::interpreter::Kind::eol);
  
  switch (token.kind) {
    case Kind::identifier:
      return token.string_value;
    case Kind::member:
      return token.string_value;
    case Kind::eof:
      return "quit";
    default:
      error->all (FC_FILE_LINE_FUNC_LINE_COL, "Expected command identifier");
      return ""; //WARNING
  }  
*/
  return "";
}

std::string Parser::get_identifier () {
  auto token = get_raw_token();
  if (token.kind == caviar::interpreter::Kind::identifier) return token.string_value;
  else error->all (FC_FILE_LINE_FUNC_LINE_COL, "Expected identifier");
  return ""; //WARNING  
}

bool Parser::compare_int () {
  return true;
}
bool Parser::compare () {
  return true;
}

bool Parser::compare_real () {
  int num_of_condition_operators = 0;
  int condition_operator = 10;
  double left_side = 0;
  double right_side = 0;
  bool operator_called = false;
  while (true) {
    auto t = get_raw_token();
    //std::cout << "found " << t.kind  << std::endl;
    //if (t.kind == caviar::interpreter::Kind::eof) break;
    //else 
    if (t.kind == caviar::interpreter::Kind::eol) break;
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
      keep_current_token();
      if (operator_called) {right_side = get_real(); break;}
      else left_side = get_real();
    } 
  }

//  std::cout <<"left: " << left_side << " ,operator: " << condition_operator << " right: " << right_side << std::endl;

  if (num_of_condition_operators == 0)
    error->all (FC_FILE_LINE_FUNC, "expected a conditional operator" ); 

  if (num_of_condition_operators > 1)
    error->all (FC_FILE_LINE_FUNC, "only one conditional operator is allowed" ); 

//  std::cout <<"left: " << left_side << " , right: " << right_side << std::endl;
///*
  if (condition_operator == -2) return (left_side <  right_side);
  if (condition_operator == -1) return (left_side <= right_side);
  if (condition_operator ==  0) return (left_side == right_side);
  if (condition_operator ==  1) return (left_side >= right_side);
  if (condition_operator ==  2) return (left_side >  right_side);

  return true;
//*/
/*
  int result_i = -5;
  int result = false;
  if (condition_operator == -2) {result = (left_side <  right_side); result_i = result;}
  if (condition_operator == -1) {result = (left_side <= right_side); result_i = result;}
  if (condition_operator ==  0) {result = (left_side == right_side); result_i = result;}
  if (condition_operator ==  1) {result = (left_side >= right_side); result_i = result;}
  if (condition_operator ==  2) {result = (left_side >  right_side); result_i = result;}
  std::cout << "result compare : " << line << " is " << result_i << std::endl;
  return result;
*/
}

std::string Parser::rest_of_line () {
  std::string st = line;
  st.erase(0, col);
  return st;
}

bool Parser::go_to_next_line () {
  while (true) {
    auto t = get_raw_token();
    if (t.kind==caviar::interpreter::Kind::eol) {
      return true;
    } else if (t.kind==caviar::interpreter::Kind::eof) {
      return false;
    }
  }
  return true;
}

int Parser::get_int () {
  return int (expression (true)); 
}

int Parser::get_positive_int () {
  int i = get_int();
  if (i<0) error->all (FC_FILE_LINE_FUNC_LINE_COL, "expected a positive integer value.");      
  return int (expression (true)); 
}

Real_t Parser::get_real () {
  return expression (true); 
}

Real_t Parser::get_positive_real () {
  Real_t i = get_real();
  if (i<0.0) error->all (FC_FILE_LINE_FUNC_LINE_COL, "expected a positive real value.");      
  return int (expression (true)); 
}

/*
// TODO: complete expression_3d before this
Vector<int> Parser::get_int_3d () {
  return int (expression_3d (true)); 
}


Vector<Real_t> Parser::get_real_3d () {
  return expression_3d (true); 
}
*/

// almost the same as 'get_literal_bool'. It should be developed so that it
// can get boolean expressions with its operators.
bool Parser::get_bool() {
  auto t = get_raw_token();
  auto s = t.string_value;
  switch (t.kind) {
    case Kind::real_number:
      return static_cast<bool>(t.real_value);
    
    case Kind::int_number:
      return static_cast<int> (t.int_value);

    case Kind::identifier:
    case Kind::string:  {
      if (s == "off" || s == "OFF" || s == "FALSE" || s == "false")
        return false;
      if (s == "on" || s == "ON" || s == "TRUE" || s == "true")
        return true;
      auto var_name = t.string_value;
      std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;
      it = object_container->dictionary.find(var_name);
      if (it == object_container->dictionary.end())   {
        error->all (FC_FILE_LINE_FUNC_LINE_COL, "unknown object name :" + var_name);      
      } else if (it->second.type == object_handler::gdst( "boolean_variable" )) {
        return object_container->boolean_variable[it->second.index];
      } else {
        std::string tmp = "boolean variable '";
        tmp += var_name;
        tmp += "' not defined";
        error->all (FC_FILE_LINE_FUNC_LINE_COL, tmp);
      } 
    }
      
    default:
      error->all (FC_FILE_LINE_FUNC_LINE_COL, "Expected boolean number literal");
  }
  return 0;//WARNING  
}

// almost the same as above
bool Parser::get_literal_bool() {
  auto t = get_raw_token();
  auto s = t.string_value;
  switch (t.kind) {
    case Kind::real_number:
      return static_cast<bool>(t.real_value);
    
    case Kind::int_number:
      return static_cast<int> (t.int_value);

    case Kind::identifier:
    case Kind::string:
      if (s == "off" || s == "OFF" || s == "FALSE" || s == "false")
        return false;
      if (s == "on" || s == "ON" || s == "TRUE" || s == "true")
        return true;
      error->all (FC_FILE_LINE_FUNC_LINE_COL, "Expected 'true' or 'false' as a bool");

      
    default:
      error->all (FC_FILE_LINE_FUNC_LINE_COL, "Expected boolean number literal");
  }
  return 0; //Warning
}

std::string Parser::get_string () {
  std::string st = "";
  Token t;
  while(true) {
    t = get_raw_token();
    if (t.kind == caviar::interpreter::Kind::eof) break; 
    else if (t.kind == caviar::interpreter::Kind::eol) break;
    else if (t.kind == caviar::interpreter::Kind::identifier) {
      auto var_name = t.string_value;
      std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;
      it = object_container->dictionary.find(var_name);
      if (it == object_container->dictionary.end())   {
        error->all (FC_FILE_LINE_FUNC_LINE_COL, "unknown object name :" + var_name);      
      } else if (it->second.type == object_handler::gdst( "string_variable" )) {
        st += object_container->string_variable[it->second.index];
      } else {
        std::string tmp = "String variable '";
        tmp += var_name;
        tmp += "' not defined";
        error->all (FC_FILE_LINE_FUNC_LINE_COL, tmp);
      } 
    } else if (t.kind == caviar::interpreter::Kind::string) {
      st += t.string_value;
    } else if (t.kind == caviar::interpreter::Kind::plus) {
      // do nothing
      //st += t.string_value;
    } else {
      error->all (FC_FILE_LINE_FUNC_LINE_COL, "expected string or string_varable");      
    }
  }
  return st; 
}

int Parser::get_literal_int () {
  auto token = get_val_token();
  if (token.kind == caviar::interpreter::Kind::int_number) return token.int_value;
  else error->all (FC_FILE_LINE_FUNC_LINE_COL, "Expected integer literal");
  return 0;//WARNING  
}

Real_t Parser::get_literal_real () {
  auto token = get_val_token();
  switch (token.kind) {
    case Kind::real_number:
      return token.real_value;
    
    case Kind::int_number:
      return token.int_value;
      
    default:
      error->all (FC_FILE_LINE_FUNC_LINE_COL, "Expected real number literal");
  }
  return 0;//WARNING  
}

std::string Parser::get_literal_string () {
  auto token = get_val_token();
  if (token.kind == caviar::interpreter::Kind::string) return token.string_value;
  else error->all (FC_FILE_LINE_FUNC_LINE_COL, "Expected string literal");
  return ""; //WARNING
}

double Parser::expression (bool get) {
  double left = term(get);
  for (;;) {
    switch (token_stream->current().kind) {
      case Kind::plus:
        left += term(true);
        break;
      case Kind::minus:
        left -= term(true);
        break;
      default:
        get_new_token = false;
        return left;
    }
  }
  return 0;//WARNING  
}

double Parser::term (bool get) {
  double left = primary(get);
  for (;;) {
    switch (token_stream->current().kind) {
      case Kind::mul:
        left *= primary(true);
        break;
      case Kind::div:
        if (auto d = primary(true)) {
          left /= d;
          break;
        }
        error->all (FC_FILE_LINE_FUNC_LINE_COL, "division by zero");
      case Kind::truediv:
        if (auto d = primary(true)) {
          left /= d;
          left = int(left);
          break;
        }
        error->all (FC_FILE_LINE_FUNC_LINE_COL, "division by zero");
      case Kind::modulus:
        if (auto d = primary(true)) {
          left = std::fmod (left, d);
          break;
        }
        error->all (FC_FILE_LINE_FUNC_LINE_COL, "division by zero");
      default:
        return left;
    }
  }
  return 0;//WARNING  
}

double Parser::primary (bool get) {
  if (get) get_raw_token ();
  switch (token_stream->current().kind) {
    case Kind::int_number: {
      int v = token_stream->current().int_value;
      get_raw_token ();
      return v;
    }
    
    case Kind::real_number: {
      double v = token_stream->current().real_value;
      get_raw_token ();
      return v;
    }

    case Kind::identifier: {
      std::string var_name = token_stream->current().string_value;
      get_raw_token ();
      std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;
      it = object_container->dictionary.find(var_name);
      if (it == object_container->dictionary.end())   {
        error->all (FC_FILE_LINE_FUNC_LINE_COL, "unknown object name :" + var_name);      
      } else if (it->second.type == object_handler::gdst( "int_variable" )) {
        return object_container->int_variable[it->second.index];
      } else if (it->second.type == object_handler::gdst( "real_variable" )) {
        return object_container->real_variable[it->second.index];
      } else {
        std::string tmp = "invalid object type '";
        tmp += var_name;
        tmp += ".";
        error->all (FC_FILE_LINE_FUNC_LINE_COL, tmp);
      } 
    }
    
    case Kind::minus:
      return -primary(true);
    
    case Kind::lp: {
      auto e = expression(true);
      if (token_stream->current().kind != Kind::rp) error->all (FC_FILE_LINE_FUNC_LINE_COL, "')' expected");
      get_raw_token ();
      return e;
    }
    
    default:
      error->all (FC_FILE_LINE_FUNC_LINE_COL, "primary expected");
  }
  return 0;//WARNING
}

// TODO: To be developed
/*
Vector<Real_t> Parser::expression_3d (bool get) {

  Vector<Real_t> left = term_3d(get);
  for (;;) {
    switch (token_stream->current().kind) {
      case Kind::plus:
        left += term_3d(true);
        break;
      case Kind::minus:
        left -= term_3d(true);
        break;
      default:
        get_new_token = false;
        return left;
    }
  }
  return 0;//WARNING  
}

Vector<Real_t> Parser::term_3d (bool get) {

  double left = primary_3d(get);
  for (;;) {
    switch (token_stream->current().kind) {
      case Kind::mul:
        left *= primary(true); // XXX Multiplying by a primary. cross_product and dot_product cannot be used here
        break;
      case Kind::div:
        if (auto d = primary(true)) { // 
          left /= d;
          break;
        }
        error->all (FC_FILE_LINE_FUNC_LINE_COL, "division by zero");
      case Kind::truediv:
        if (auto d = primary_3d(true)) {//XXX Division by an primary not an primary_3d
          left /= d;
          left = int(left);
          break;
        }
        error->all (FC_FILE_LINE_FUNC_LINE_COL, "division by zero");
      case Kind::modulus:
        if (auto d = primary_3d(true)) {
          left = std::fmod (left, d);
          break;
        }
        error->all (FC_FILE_LINE_FUNC_LINE_COL, "division by zero");
      default:
        return left;
    }
  }
  return 0;//WARNING  

}

Vector<Real_t> Parser::primary_3d (bool get) {

  if (get) get_raw_token ();
  switch (token_stream->current().kind) {
    case Kind::int_number: {
      int v = token_stream->current().int_value;
      get_raw_token ();
      return v;
    }
    
    case Kind::real_number: {
      double v = token_stream->current().real_value;
      get_raw_token ();
      return v;
    }

    case Kind::identifier: {
      std::string var_name = token_stream->current().string_value;
      get_raw_token ();
      std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it;
      it = object_container->dictionary.find(var_name);
      if (it == object_container->dictionary.end())   {
        error->all (FC_FILE_LINE_FUNC_LINE_COL, "unknown object name :" + var_name);      
      } else if (it->second.type == object_handler::gdst( "int_variable" )) {
        return object_container->int_variable[it->second.index];
      } else if (it->second.type == object_handler::gdst( "real_variable" )) {
        return object_container->real_variable[it->second.index];
      } else {
        std::string tmp = "invalid object type '";
        tmp += var_name;
        tmp += ".";
        error->all (FC_FILE_LINE_FUNC_LINE_COL, tmp);
      } 
    }
    
    case Kind::minus:
      return -primary(true);
    
    case Kind::lp: {
      auto e = expression(true);
      if (token_stream->current().kind != Kind::rp) error->all (FC_FILE_LINE_FUNC_LINE_COL, "')' expected");
      get_raw_token ();
      return e;
    }
    
    default:
      error->all (FC_FILE_LINE_FUNC_LINE_COL, "primary expected");
  }
  return 0;//WARNING

}
*/
} //interpreter
} // namespace caviar

