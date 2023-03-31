
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
//#include "caviar/objects/atom_data.h"
//#include <ctime>

namespace caviar {

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
    if (string_cmp(t,"output_all_acc")) {
      write();
    }/* else  if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } */
  }
  return in_file;
}

void Force_field::initialize(){}
void Force_field::write(){
  /*std::ofstream ofs ("o_acc");
  const auto &pos = atom_data -> owned.position;  
  const auto &acc = atom_data -> owned.acceleration;  
  for (unsigned int i=0;i<pos.size();++i) {
    ofs << i << " " << acc[i].x << "\t" << acc[i].y << "\t" << acc[i].z << "\n" ;
  }*/
}
void Force_field::write(int64_t){} // current time_step
void Force_field::write(double){} // current time
void Force_field::write(int64_t, double){} //time_step and time
void Force_field::start_new_files(){} //add_time_to_previous
void Force_field::start_new_files(std::string &){} //add_time_to_previous

} //Force_field

} // namespace caviar

