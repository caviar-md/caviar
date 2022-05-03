
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

#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/interpreter/object_handler/preprocessors_new.h"

#include <string>

namespace caviar {
namespace objects {

bool Atom_data::add_xyz_data_file (caviar::interpreter::Parser *parser) {
  output->info("Basic::add_xyz_data_file ");

  std::string xyz_file_name = "";
  bool last_frame = false;
  bool in_file = true;
  bool read_velocity = false;
  bool replace_data = false; // used for resuming simulations when there's bonds and angles
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"last_frame")) {
      last_frame = true;
    } else if (string_cmp(t,"replace_data")) {
      replace_data = true;
    } else if (string_cmp(t,"file_name")) {
      const auto token = parser->get_val_token();
      const auto file_name = token.string_value;
//      GET_A_STRING(xyz_file_name,"","")
      xyz_file_name = file_name;
    } else if (string_cmp(t,"read_velocity")) {
      read_velocity = true;
    } else FC_ERR_UNDEFINED_VAR(t)

  }

  int start_line = 1;

  if (last_frame) {
    int i = 1;
    int num_xyz_frames = 0;
//    Parser *pf (xyz_file_name);
    caviar::interpreter::Parser pf (fptr, xyz_file_name);
    auto t = pf.get_val_token();
    while (t.kind != caviar::interpreter::Kind::eof) {
      if (t.kind == caviar::interpreter::Kind::identifier) {
        auto ts = t.string_value;
        if (   string_cmp(ts,"atom") || string_cmp(ts,"atoms")
            || string_cmp(ts,"ATOM") || string_cmp(ts,"ATOMS") 
        ) {
          ++num_xyz_frames;
          if (i==1) error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown xyz format");
          start_line = i-1;
        }
      }
      pf.end_of_line();
      t = pf.get_val_token();
      ++i;
    }
    std::string st1 = "Num. of total xyz frames found : " + std::to_string(num_xyz_frames);    
    std::string st2 = "Last xyz frame found at line number : " + std::to_string(i);
    output->info(st1);
    output->info(st2);
  }


  caviar::interpreter::Parser pf (fptr, xyz_file_name);
  for (int i = 1; i <start_line; ++i) {
    pf.end_of_line();
  }

/*
  int num_atoms = pf.get_literal_int();
  std::cout << "num_atoms : " << num_atoms << std::endl;
      pf.end_of_line();
  auto st = pf.get_val_token();
  std::cout << "st : " << st.string_value << std::endl;
      pf.end_of_line();
*/

  int num_atoms = pf.get_literal_int();

  pf.end_of_line();

  std::string st = " Number of atoms found at the inputted xyz file : " + std::to_string(num_atoms) ;

  if (num_atoms==0) error->all(FC_FILE_LINE_FUNC, "no atom found in the xyz file. Something is wrong.");

  output->info(st);
  pf.get_val_token();
  pf.end_of_line();

  if (replace_data) {
    if (num_atoms != (int) owned.position.size()) error->all (FC_FILE_LINE_FUNC, "XYZ file is not compatible: Different number of existing atoms in atomdata and xyz file atoms.");
  }
  
  for (int i = 0; i <num_atoms; ++i) {
    auto type = pf.get_literal_int();

    Vector<Real_t> pos, vel{0.0, 0.0, 0.0};;

    pos.x = pf.get_literal_real();
    pos.y = pf.get_literal_real();
    pos.z = pf.get_literal_real();

    if (read_velocity) {
      vel.x = pf.get_literal_real();
      vel.y = pf.get_literal_real();
      vel.z = pf.get_literal_real();
    }

    pf.end_of_line();

    
    
    if (replace_data) {
        if (type != (int) owned.type[i]) error->all (FC_FILE_LINE_FUNC, "XYZ file is not compatible: Different atom type order exists.");
        owned.position[i] = pos;
        owned.velocity[i] = vel;
    } else {
        auto id = get_global_id();
        add_atom (id, type, pos, vel);
    }

  }

  return in_file;
}


} //objects
} // namespace caviar

