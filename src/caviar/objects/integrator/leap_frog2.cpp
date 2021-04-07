
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
#include "caviar/objects/integrator/leap_frog2.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/force_field.h"

namespace caviar {
namespace objects {
namespace integrator {

Leap_frog2::Leap_frog2 (CAVIAR *fptr) : Integrator{fptr} {
  FC_OBJECT_INITIALIZE_INFO
  integrator_type = 2;
}

Leap_frog2::~Leap_frog2 (){}

bool Leap_frog2::read (caviar::interpreter::Parser *parser) {
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

void Leap_frog2::verify_settings (){
  FC_NULLPTR_CHECK(atom_data)
}

// Rapaport 'Art of Molecular Dynamics'. The first leap-frog scheme.
// Is this implementation OK? XXX 
void Leap_frog2::step_part_I () {

  FC_OBJECT_VERIFY_SETTINGS

  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;

  const auto psize = vel.size();

#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 
    vel [i] += dt * acc[i];             // v (t + dt/2) = v  (t − dt/2) + dt * a  (t) // XXX
  }
}

void Leap_frog2::step_part_II () {

  FC_OBJECT_VERIFY_SETTINGS

  auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;

  const auto psize = pos.size();

#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 
    pos [i] +=  dt * vel [i];         // r (t + dt) = r (t) + dt * v (t + dt/2)
  }
}

void Leap_frog2::step_part_III () {
  // NOTHING
}

} //integrator
} //objects

} // namespace caviar


