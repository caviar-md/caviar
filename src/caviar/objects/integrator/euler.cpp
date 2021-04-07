
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

#include "caviar/objects/integrator/euler.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/force_field.h"

namespace caviar {
namespace objects {
namespace integrator {

Euler::Euler (CAVIAR *fptr) : Integrator{fptr} {
  FC_OBJECT_INITIALIZE_INFO
  integrator_type = 4;
}

Euler::~Euler (){}

bool Euler::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"dt")) {
      GET_OR_CHOOSE_A_REAL(dt,"","")
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  return in_file;
}

void Euler::verify_settings (){
  FC_NULLPTR_CHECK(atom_data)
  atom_data->record_owned_acceleration_old = true;
}

void Euler::step_part_I () {
  // NOTHING
}

void Euler::step_part_II () {

  FC_OBJECT_VERIFY_SETTINGS

  auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;

  const auto psize = pos.size();
  
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 
    pos [i] += vel [i] * dt + 0.5 * acc [i] * dt * dt; // r (t+dt) = r(t) + v(t)*dt + 1/2 * a(t) * dt^2
    vel [i] += acc [i] * dt;                           // v (t+dt) = v(t) + a(t)*dt;

  }


}


void Euler::step_part_III () {
  // NOTHING
}

} //integrator
} //objects

} // namespace caviar


