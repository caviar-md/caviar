
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

#include "caviar/interpreter/input/lexer.h"
#include "caviar/interpreter/communicator.h"
#include "caviar/interpreter/error.h"

//#define DEBUG_ME

CAVIAR_NAMESPACE_OPEN
namespace interpreter {
static constexpr size_t max_buffer_size = 1024;

Token_stream::Token_stream (CAVIAR *fptr) : Pointers{fptr},
    input_stream{&std::cin}, stream_is_file{false}, get_new_line{true},
    multiline_comment{false} {
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: LEXER constructor t1" << std::endl;
#endif
}

Token_stream::Token_stream (CAVIAR *fptr, const std::string &file) : Pointers{fptr},
    input_stream{new std::ifstream {file}}, stream_is_file{true},
    get_new_line{true}, multiline_comment{false}   {
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: LEXER constructor t2" << std::endl;
#endif
  if (input_stream->fail()) {
    std::string tmp = "Can not open file: '";
    tmp += file;
    tmp += '\'';
    error->all (FC_FILE_LINE_FUNC_LINE_COL, tmp);
  }
}

Token_stream::Token_stream (CAVIAR *fptr, std::istringstream &iss) : Pointers{fptr}, input_stream{&iss}, stream_is_file{false}, get_new_line{true} {
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: LEXER constructor t3" << std::endl;
#endif
}

Token_stream::~Token_stream () {
#ifdef DEBUG_ME
  std::cout<<"DEBUG_ME: LEXER destructor t1" << std::endl;
#endif
  if (stream_is_file) delete input_stream;
}

void Token_stream::getline () {
  col = 0;
#if defined(CAVIAR_WITH_MPI) || defined(CAVIAR_WITH_DEALII_MPI)
  if (comm->me == 0) {
    std::getline (*input_stream, line);
    end_of_file = input_stream->eof();
    line += '\n';
  } else
    line = "";
  comm->broadcast (end_of_file);
  comm->broadcast (line);
#else

  std::getline (*input_stream, line);
  //std::cout << "LLX: " << line << std::endl;
  end_of_file = input_stream->eof();
  //if (end_of_file)
    //std::cout << "Wow : EOF found." << std::endl;

  line += '\n';
#endif
  line_stream = std::istringstream (line);
  get_new_line = false;
}

char Token_stream::get_nonblank_char () {
  char ch;
  do {
    ch = get_char ();
  } while (isblank (ch));
  return ch;
}

//  modulus='%', assign='=',
//  lp='(', rp=')', squote='\'', dquote='"', comment='#', backslash='\\', eol='\n',
//  lsb='[', rsb=']', lcb='{', rcb='}', larger='>',  smaller='<',
//  colon=':', comma=',',


Token Token_stream::get () {

  ct = Kind::unknown;

  if (get_new_line) getline();

  if (end_of_file) 
    return ct={Kind::eof};

  char ch = get_nonblank_char();
  
  while (ch == '\\') {
    if (get_char () == '\n') {
      getline();
      ch = get_nonblank_char();
    } else error->all (FC_FILE_LINE_FUNC_LINE_COL, "unexpected character after line continuation character");
  }
  
  // XXX: For the futue: Multiline comment has meaning when there is multiline 
  // commands. But in one command per line case, there may be exceptions.
  // For example when one wants to assign 'int i = 3' with a comment,
  // he or she has many options as the two cases below. The question is that 
  // is 'case 2' a valid caviar command?
  //
  // case 1: 
  //   int i #* this is a counter variable *# = 3
  //
  // case 2:
  //   int i #* this is a 
  //            counter variable *# = 3
  //
  // 
  /*
  if (multiline_comment) {
    switch (ch) {
      case '\n':
        get_new_line = true;
        return ct = {Kind::comment};
      
      case '*': 
        ch = get_char ();
        if (ch == '#') {
          multiline_comment = false;
        } else {
          putback_char (ch);
        }
        return ct = {Kind::comment};

      default:
        return ct = {Kind::comment};
    }
  }   
  */

  switch (ch) {
    case '%':
      return ct = {Kind::modulus};

    case '(':
      return ct = {Kind::lp};

    case ')':
      return ct = {Kind::rp};
 
    case '[':
      return ct = {Kind::lsb};

    case ']':
      return ct = {Kind::rsb};

    case '{':
      return ct = {Kind::lcb};

    case '}':
      return ct = {Kind::rcb};

    case ':': 
      return ct = {Kind::colon};

    case ',':
      return ct = {Kind::comma};

    case '+':
      ch = get_char ();
      if (ch == '=') {
        return ct = {Kind::asplus};
      }
      else {
        putback_char (ch);
        return ct = {Kind::plus};
      }

    case '-':
      ch = get_char ();
      if (ch == '=') {
        return ct = {Kind::asminus};
      }
      else {
        putback_char (ch);
        return ct = {Kind::minus};
      }

    case '*':
      ch = get_char ();
      if (ch == '=') {
        return ct = {Kind::asmul};
      }
      else {
        putback_char (ch);
        return ct = {Kind::mul};
      }

    case '=':
      ch = get_char ();
      if (ch == '=') {
        return ct = {Kind::equal};
      }
      else {
        putback_char (ch);
        return ct = {Kind::assign};
      }

    case '>':
      ch = get_char ();
      if (ch == '=') {
        return ct = {Kind::eqlarger};
      }
      else {
        putback_char (ch);
        return ct = {Kind::larger};
      }    

    case '<':
      ch = get_char ();
      if (ch == '=') {
        return ct = {Kind::eqsmaller};
      }
      else {
        putback_char (ch);
        return ct = {Kind::smaller};
      }

    case '/':
      ch = get_char ();
      if (ch == '/') {
        return ct = {Kind::truediv};
      } else if (ch == '='){
        return ct = {Kind::asdiv};
      } else {
        putback_char (ch);
        return ct = {Kind::div};
      }
    
/*
    case '#':
      ch = get_char ();
      if (ch == '*') {  // multiline comment start point
        multiline_comment = true;
        return ct = {Kind::comment};
      } else {          // single line comment
        putback_char (ch);
        get_new_line = true;
        return ct = {Kind::eol};
      }
*/
    case '#':
    case '\n':
      get_new_line = true;
      return ct = {Kind::eol};

    
    case '.':
    case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
      return get_number_literal (ch);
    
    case '\'':
    case '"':
      return get_string_literal (ch);
    
    // get identifier or member object
    default:
      if (isalpha(ch) || ch == '_') {
        ct.string_value = ch;
        int member_index = -1;
        while ((ch = get_char ()) && (isalnum(ch) || ch == '_' || ch == '.')) {
          if (ch == '.') {
            ++member_index;
            ct.member.push_back("");
            continue;
          }
          if (member_index == -1)
            ct.string_value += ch; // append ch to end of string_value
          else
            ct.member[member_index] += ch; 
        }

        putback_char (ch);
        if (member_index == -1)
          ct.kind = Kind::identifier;
        else
          ct.kind = Kind::member;
        return ct;
      }
      error->all (FC_FILE_LINE_FUNC_LINE_COL, "Invalid token");
  }
  return ct;//WARNING
}

Token & Token_stream::current () {
  return ct;
}

Token Token_stream::get_number_literal (char first_char) {
  bool num_is_real = (first_char == '.');
  std::string tmp;
  tmp += first_char;
  char ch;
  
  while (!num_is_real) {
    ch = get_char ();
    tmp += ch;
    if (isdigit(ch)) continue;
    else if (ch == '.') {
      num_is_real = true;
      break;
    } else if (ch == 'e' || ch == 'E') {
      ch = get_char ();
      tmp += ch;
      if (isdigit(ch)) {
        num_is_real = true;
        break;
      } else if (ch == '+' || ch == '-') {
        ch = get_char ();
        tmp += ch;
        if (isdigit(ch)) {
          num_is_real = true;
          break;
        } else error->all (FC_FILE_LINE_FUNC_LINE_COL, "Invalid token");
      } else error->all (FC_FILE_LINE_FUNC_LINE_COL, "Invalid token");
    } else break;
  }
  
  while (tmp.length()) {
    putback_char (tmp.back());
    tmp.pop_back();
  }
  
  ct.kind = num_is_real ? Kind::real_number : Kind::int_number;
  if (num_is_real) line_stream >> ct.real_value;
  else line_stream >> ct.int_value;
  return ct;
}

Token Token_stream::get_string_literal (char delim) {
  ct.kind = Kind::string;
  ct.string_value.clear();
  char ch, ch2, ch3;
  bool triple = false;
  constexpr auto eof = std::char_traits<char>::eof();
  
  if ((ch = get_char ()) == delim) {
    if ((ch2 = get_char ()) == delim) triple = true;
    else {
      putback_char (ch2);
      putback_char (ch);
    }
  } else putback_char (ch);
  
  while (line_stream.get(ch)) {
    if (ch == '\\') {
      switch (get_char ()) {
        case '\'':
          ct.string_value += '\'';
          break;
        case '"':
          ct.string_value += '"';
          break;
        case '\n':
          getline();
          break;
        case 'a':
          ct.string_value += '\a';
          break;
        case 'b':
          ct.string_value += '\b';
          break;
        case 'f':
          ct.string_value += '\f';
          break;
        case 'n':
          ct.string_value += '\n';
          break;
        case 'r':
          ct.string_value += '\r';
          break;
        case 't':
          ct.string_value += '\t';
          break;
        case 'v':
          ct.string_value += '\v';
          break;
        case '\\':
        default:
          ct.string_value += '\\';
      }
    } else if (ch == delim) {
      if (triple) {
        if ((ch2 = get_char ()) == delim) {
          if ((ch3 = get_char ()) == delim) break;
          else {
            if (ch3 == '\n') getline();
            ct.string_value += ch;
            ct.string_value += ch2;
            ct.string_value += ch3;
          }
        } else {
          if (ch2 == '\n') getline();
          ct.string_value += ch;
          ct.string_value += ch2;
        }
      } else break;
    } else if (ch == '\n') {
      if (triple) {
        ct.string_value += '\n';
        getline();
      } else error->all (FC_FILE_LINE_FUNC_LINE_COL, "EOL while scanning string literal");
    } else if (ch == eof) error->all (FC_FILE_LINE_FUNC_LINE_COL, "EOF while scanning triple-quoted string literal");
    else ct.string_value += ch;
  }
  
  return ct;
}
} //interpreter
CAVIAR_NAMESPACE_CLOSE

