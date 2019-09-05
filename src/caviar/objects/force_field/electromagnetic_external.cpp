
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

#include "caviar/objects/force_field/electromagnetic_external.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"

#include <cmath>
#include <iomanip>

namespace caviar {
namespace objects {
namespace force_field {

Electromagnetic_external::Electromagnetic_external (CAVIAR *fptr) : Force_field{fptr},
    amplitude_E{1.0}, amplitude_B{1.0}, direction_E{Vector<double>{0,0,0}},
    direction_B{Vector<double>{0,0,0}}
{
  FC_OBJECT_INITIALIZE_INFO
}

bool Electromagnetic_external::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"amplitude_E")) {
      GET_OR_CHOOSE_A_REAL(amplitude_E,"","")
    } else if (string_cmp(t,"direction_E")) {
      GET_OR_CHOOSE_A_REAL_3D_VECTOR(direction_E, "", "");
      auto d_sq = direction_E*direction_E;
      auto d_norm = std::sqrt(d_sq);
      direction_E = direction_E/d_norm;
    } else if (string_cmp(t,"amplitude_B")) {
      GET_OR_CHOOSE_A_REAL(amplitude_B,"","")
    } else if (string_cmp(t,"direction_B")) {
      GET_OR_CHOOSE_A_REAL_3D_VECTOR(direction_B, "", "");
      auto d_sq = direction_B*direction_B;
      auto d_norm = std::sqrt(d_sq);
      direction_B = direction_B/d_norm;
    } else if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  
  return in_file;
}

void Electromagnetic_external::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)
}


void Electromagnetic_external::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS


  const auto &pos = atom_data -> owned.position;  
  const auto &vel = atom_data -> owned.position;

  for (unsigned int i=0;i<pos.size();++i) {
    const auto type_i = atom_data -> owned.type [i] ;
    const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];
    const auto charge_i = atom_data -> owned.charge [ type_i ];


    const auto a = charge_i * (amplitude_E * direction_E + amplitude_B * cross_product(vel[i], direction_B)) * mass_inv_i;
    atom_data -> owned.acceleration [i] += a;
    
  }  

}

} //force_field
} //objects
} // namespace caviar

