
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

#include "caviar/objects/writer/force_field.h"
#include "caviar/utility/interpreter_io_headers.h"
//#include <ctime>

namespace caviar {
namespace objects {
namespace writer {

Force_field::Force_field (CAVIAR *fptr) : Writer{fptr}
{
  FC_OBJECT_INITIALIZE_INFO
}

Force_field::~Force_field () {}

bool Force_field::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;

  }
  return in_file;
}

void Force_field::initialize(){}
void Force_field::write(){}
void Force_field::write(int){} // current time_step
void Force_field::write(double){} // current time
void Force_field::write(int, double){} //time_step and time
void Force_field::start_new_files(){} //add_time_to_previous
void Force_field::start_new_files(std::string &){} //add_time_to_previous

} //Force_field
} //objects
} // namespace caviar

