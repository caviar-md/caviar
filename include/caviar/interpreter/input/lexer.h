
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

#ifndef CAVIAR_INTERPRETER_LEXER_H
#define CAVIAR_INTERPRETER_LEXER_H

#include <string>
#include <istream>
#include <sstream>
#include <vector>

#include "caviar/utility/pointers.h"

namespace caviar {
namespace interpreter {
// enum class's underlying type is 'char'
// 'enum class' is preferred over 'enum' since the names are local.
/*
enum class Kind : char {
  identifier, string, int_number, real_number, eof,
  plus='+', minus='-', mul='*', div='/', truediv, modulus='%', assign='=', asplus, asminus, asmul, asdiv,
  lp='(', rp=')', squote='\'', dquote='"', comment='#', backslash='\\',eol='\n',
  lsb='[', rsb=']', lcb='{', rcb='}', larger='>', eqlarger, smaller='<',
  eqsmaller, equal, colon=':', comma=','
};
*/

/**
 * 'Kind' being a char is not neccesary. It was used like above with some use of
 * static_cast<Kind> (SOME_CHAR) in lexer functions. But since it arises some 
 * problem in their underlaying char value, i.e., because of the sequence of the
 * char is not neccesarily obeyed there would be some Kind which are the same,
 * we changed it to this. For reference, we keep the above scheme in the comments.
 */
enum class Kind : char {
  eqsmaller, equal,
  identifier, string, int_number, real_number, eof,
  asplus, asminus, asmul, asdiv,
  truediv, eqlarger,
  plus, minus, mul, div,  modulus, assign,
  lp, rp, squote, dquote, comment, backslash, eol,
  lsb, rsb, lcb, rcb, larger,  smaller,
  colon, comma, member, unknown
};

/**
 * Tokens structure and their variables
*/
struct Token {
  Token () {}
  Token (Kind kind) : kind{kind} {}
  Kind kind;
  // making this a union will complicate destructor and assignment operator b/c of std::string
  std::string string_value;
  std::vector<std::string> member;
  int int_value;
  Real_t real_value;
};

/**
 * a class that gets a string stream and tokenize it for the parser.
 */
class Token_stream : public Pointers {
public:
  // set's input_stream to std::cin
  Token_stream (class CAVIAR *);
  // set's input_stream to std::ifstream of a file
  Token_stream (class CAVIAR *, const std::string &);
  // set's input_stream to std::stringstream. 
  Token_stream (class CAVIAR *, std::istringstream &);
  ~Token_stream ();
  Token get ();

  Token & current ();
  void set_input (std::istream &);
  std::string line;
  unsigned int col;
public:
  void getline ();
  inline char get_char () {++col; return line_stream.get ();}
  inline void putback_char (char ch) {--col; line_stream.putback (ch);}
  char get_nonblank_char ();
  Token get_number_literal (char);
  Token get_string_literal (char);
  
  std::istream *input_stream;
  std::istringstream line_stream;
  Token ct;
  bool stream_is_file, get_new_line, end_of_file;
  bool multiline_comment;
};

} //interpreter
} // namespace caviar

#endif
