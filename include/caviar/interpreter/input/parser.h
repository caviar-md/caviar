
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

#ifndef CAVIAR_INTERPRETER_PARSER_H
#define CAVIAR_INTERPRETER_PARSER_H

#include "caviar/utility/pointers.h"
#include "caviar/utility/vector.h"

namespace caviar {
namespace interpreter {
/**
 * Parser, is an important part of script reading by any code that does so.
 * 'Input' class creates a 'Parser' out of a file or an input stream, then gets
 * the input from the parser.
 * Every 'Parser' creates a 'Token_stream' as it initializes. 'Token_stream'
 * tokenize the input stream (in most cases, an input file) line by line and gives 
 * it to the parser. Parser organizes these tokens into 'commands', 'mathematical
 * expressions' and 'variables'. Running basic commands is what is 'Input' class
 * does. Input class also gives unknown commands to 'Object_creator', 'Object_handler'
 * and 'Object_container' classes to check if it is their Command. 
 */
class Parser : public Pointers {
public:
  Parser (class CAVIAR *);
  Parser (class CAVIAR *, const std::string &);
  Parser (class CAVIAR *, std::istringstream &);

  ~Parser ();

  /**
   * Get a raw token and check if it is the end of line.
   */
  bool end_of_line (); 

  /**
   * Get a raw token and check if it is assignment, i.e., "="
   */
  bool assignment (); 

  /**  
   * Affects getting a new token in 'get_raw_token()' function.
   */
  void keep_current_token () {get_new_token = false;}

  /**
   * Just gets a raw token! As it says :)
   */
  struct Token get_raw_token (); 

  /**
   * Gets a token. In case of '+' or '-' gets the next token. If it's not
   * a number, returns an error.
   */
  struct Token get_val_token (); 
  
  /**
   * Write a raw string or a string variable into 'std::ostream'.
   */
  std::ostream & write_to_stream (std::ostream &); 

  /**
   * if a string is in '"', it is returned as a Kind::string. If not, it will
   * be returned as an Kind::identifier. Note that every identifier has a 
   * string_value.
   * Gets a identifier token as a std::string. If not, returns an error.
   */
  std::string get_identifier (); 

  /**
   * Not much different from above, only used once in 'input.cpp'.
   * Needs some illustrations.
   */
  std::string get_command_identifier (); 

  /**
   * returns an int value out of an expression with mathematical operators and 
   * numeric variables.
   */
  int get_int ();
  int get_positive_int();

  /**
   * returns a real value out of an expression with mathematical operators and 
   * numeric variables.
   */
  Real_t get_real (); 
  Real_t get_positive_real (); 

  /**
   * parses the conditional line to 'eol' and gives the result.
   */
  bool compare_real ();
  bool compare_int ();
  bool compare ();

  /**  
   * gets a bool value of an expression . TODO: complete it
   */
  bool get_bool ();

  /**
   * gets a bool, such as 'true', 'false', 'on', 'off'
   */
  bool get_literal_bool();

  //Vector<Real_t> get_real3d ();
  //Vector<int> get_int3d ();

  /**
   * returns a string. It can contains string variables. 
   */
  std::string get_string ();

  /**  
   * Gets a int value in literal numbers. No variable.
   */
  int get_literal_int (); 

  /**
   * gets a real value in literal numbers. No variable!
   */
  Real_t get_literal_real ();

  /**
   * gets a string value. It should be in qouts, i.e. " something "
   */
  std::string get_literal_string ();

  /**
   * give the rest of line, line after where it is parsed to 'eol'
   */
  std::string rest_of_line ();

  bool go_to_next_line ();

public:
  class Token_stream *token_stream;
  bool get_new_token;

  /**
   * Every 'expression' contains some 'term' which can have some 'primary'.
   * every primary can be an expression itself, using '(' .
   * int expressions are also calculated by these functions.
   */
  double expression (bool);
  double term (bool);
  double primary (bool);

  /**
   * CAVIAR::Vector counterpart of the functions above
   */
  //Vector<Real_t> expression_3d (bool);
  //Vector<Real_t> term_3d (bool);
  //Vector<Real_t> primary_3d (bool);
  
public:

  /**
   * current 'line' and 'column' of the token_stream.
   */
  std::string &line;
  unsigned int &col;
};
} //interpreter
} // namespace caviar

#endif
