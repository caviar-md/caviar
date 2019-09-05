
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

#include "caviar/objects/force_field/gravity_external.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include <cmath>
#include <iomanip>

namespace caviar {
namespace objects {
namespace force_field {

Gravity_external::Gravity_external (CAVIAR *fptr) : Force_field{fptr},
 amplitude{1.0}, direction{Vector<double>{0,0,0}}
 {
  FC_OBJECT_INITIALIZE_INFO
}

bool Gravity_external::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"amplitude")) {
      GET_OR_CHOOSE_A_REAL(amplitude,"","")
    } else if (string_cmp(t,"direction")) {
      GET_OR_CHOOSE_A_REAL_3D_VECTOR(direction, "", "");
      auto d_sq = direction*direction;
      auto d_norm = std::sqrt(d_sq);
      direction = direction/d_norm;
    } else if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  
  return in_file;
}

void Gravity_external::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)
}


void Gravity_external::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS


  const auto &pos = atom_data -> owned.position;  

  for (unsigned int i=0;i<pos.size();++i) {
    //const auto type_i = atom_data -> owned.type [i] ;
    //const auto mass_i = atom_data -> owned.mass [ type_i ];

    //const auto force = amplitude * direction * mass_i;
    //atom_data -> owned.acceleration [i] += force / mass_i;

    const auto a = amplitude * direction;
    atom_data -> owned.acceleration [i] += a;
    
  }  

}

} //force_field
} //objects
} // namespace caviar

