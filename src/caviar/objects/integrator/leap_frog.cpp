
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
#include "caviar/objects/integrator/leap_frog.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/force_field.h"

namespace caviar {
namespace objects {
namespace integrator {

Leap_frog::Leap_frog (CAVIAR *fptr) : Integrator{fptr} {
  FC_OBJECT_INITIALIZE_INFO
  integrator_type = 2;
}

Leap_frog::~Leap_frog (){}

bool Leap_frog::read (caviar::interpreter::Parser *parser) {
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
    } else  FC_ERR_UNDEFINED_VAR(t)
  }
  return in_file;
}

void Leap_frog::verify_settings (){
  FC_NULLPTR_CHECK(atom_data)
}

void Leap_frog::step_part_I () {
  FC_OBJECT_VERIFY_SETTINGS

  auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;

  const auto psize = pos.size();

#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 
    vel [i] += 0.5 * dt * acc[i]; // v (t + dt/2) = v (t) + (dt/2) a (t)
  }


}

void Leap_frog::step_part_II () {

  auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;
  //auto &acc = atom_data -> owned.acceleration;

  const auto psize = pos.size();

#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 
    pos [i] += dt * vel[i];       // r (t + dt) = r (t) + dt * v (t + dt/2)
  }


}

void Leap_frog::step_part_III () {
  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;
  
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<vel.size(); i++) {
    vel [i] += 0.5 * dt * acc [i]; // v (t + dt) = v (t + dt/2) + (dt/2) a (t + dt)
  }

}

} //integrator
} //objects

} // namespace caviar


