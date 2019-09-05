
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

#include "caviar/objects/integrator/verlet.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"
#include "caviar/objects/force_field.h"

namespace caviar {
namespace objects {
namespace integrator {

Verlet::Verlet (CAVIAR *fptr) : Integrator{fptr} {
  FC_OBJECT_INITIALIZE_INFO
  integrator_type = 4;
}

Verlet::~Verlet (){}

bool Verlet::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"set_domain") || string_cmp(t,"domain")) {
      FIND_OBJECT_BY_NAME(domain,it)
      domain = object_container->domain[it->second.index];
    } else if (string_cmp(t,"dt")) {
      GET_OR_CHOOSE_A_REAL(dt,"","")
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  return in_file;
}

void Verlet::verify_settings (){
  FC_NULLPTR_CHECK(atom_data)
  FC_NULLPTR_CHECK(domain)
  atom_data->record_owned_position_old = true;
}

// This is implemented according to
// 'Computational Physics - Molecular Dynamics Simulations-
// E. Carlon, M. Laleman and S. Nomidis – Academic year 2015/2016'
// and also,
// 'An overview of integration schemes for molecular dynamics simulations-
// Ulf D. Schiller - 5th March 2008'
//

void Verlet::step_part_I () {
  // nothing to do here.
}

void Verlet::step_part_II () {

  FC_OBJECT_VERIFY_SETTINGS

  auto &pos = atom_data -> owned.position;
  auto &pos_old = atom_data -> owned.position_old;
  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;

  const auto psize = pos.size();

  const auto two_dt_inv = 1.0/(2.0*dt);

  for (unsigned int i=0; i<psize; i++) { 
    // r (t + dt) = 2*r(t) − r (t − dt) + a(t) * dt^2 
    pos [i] = (2.0*pos[i]) - pos_old[i]  + acc [i] * dt * dt; 

    // v(t) = (r (t + dt) - r (t - dt) )/ (2*dt) 
    // XXX note that this is 'v(t)' not 'v(t+dt)'
    vel [i] = domain->fix_distance(pos[i] - pos_old[i]) * two_dt_inv; 

  }

}


void Verlet::step_part_III () {
  // NOTHING
}

} //integrator
} //objects

} // namespace caviar


