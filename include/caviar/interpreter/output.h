
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

#ifndef CAVIAR_INTERPRETER_OUTPUT_H
#define CAVIAR_INTERPRETER_OUTPUT_H

#include "caviar/utility/pointers.h"

namespace caviar {
namespace interpreter {

/**
 * This class does all of the output massages.
 * 
 * 
 */
class Output : public Pointers {
public:
  Output (class CAVIAR *);

  void comment (const std::string &, const bool endline=true); 
  void comment (const char *, const bool endline=true); 

  void info (const std::string &, const int level=0, const bool endline=true);
  void info (const char *, const int level=0, const bool endline=true);

  void info_create (const std::string &, const int level=1, const bool endline=true);
  void info_create (const char *, const int level=1, const bool endline=true);

  void info_read (const std::string &, const int level=1, const bool endline=true);
  void info_read (const char *, const int level=1, const bool endline=true);

  void warning (const std::string &, const int level=0, const bool endline=true);
  void warning (const char *, const int level=0, const bool endline=true);

  bool read (class caviar::interpreter::Parser *);

  bool output_info[5], output_warning[5];

};


} //interpreter
} // namespace caviar

#endif
