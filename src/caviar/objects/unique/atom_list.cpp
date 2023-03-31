
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

#include "caviar/objects/unique/atom_list.h"
#include "caviar/objects/unique/atom.h"
#include "caviar/objects/unique/atom_group.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/interpreter/object_handler/preprocessors_new.h"

namespace caviar {

namespace unique {

Atom_list::Atom_list (CAVIAR *fptr) : Unique{fptr} {
  FC_OBJECT_INITIALIZE_INFO
}   

   
Atom_list::~Atom_list () {}


void Atom_list::verify_settings () {
  
}


Atom_list::Atom_list (const Atom_list & a) : Unique{a} {

}


bool Atom_list::read ( caviar::interpreter::Parser * parser) {
  FC_OBJECT_READ_INFO

  while(true) {
    FC_IF_RAW_TOKEN_EOF_EOL
    if (string_cmp(ts,"add_atom")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,atom)
      auto a =  dynamic_cast<unique::Atom *>(object_container->unique[it->second.index]);

      atoms.push_back(a);
      continue;
    } else if (string_cmp(ts,"clear")) {
      atoms.clear(); continue;
    } else FC_ERR_UNDEFINED_VAR(ts)    
  }

  return true;
}

void Atom_list::add_atom(const unique::Atom &) {

}

} //unique


} // namespace caviar

