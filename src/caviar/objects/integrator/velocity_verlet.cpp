
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

#include "caviar/objects/integrator/velocity_verlet.h"
//#include "caviar/utility/python_utils_def.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/force_field.h"

namespace caviar {
namespace objects {
namespace integrator {

Velocity_verlet::Velocity_verlet (CAVIAR *fptr) : Integrator{fptr} {
  FC_OBJECT_INITIALIZE_INFO
  integrator_type = 1;
  type = 1;
}

Velocity_verlet::~Velocity_verlet (){}

bool Velocity_verlet::read (caviar::interpreter::Parser *parser) {
  /*
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
    } else if (string_cmp(t,"type")) {
      GET_OR_CHOOSE_A_INT(type,"","")
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  return in_file;
  */
}

void Velocity_verlet::verify_settings (){
  FC_NULLPTR_CHECK(atom_data)

  switch (type) {


  case 1:
  break;

  case 2:{
    atom_data->record_owned_acceleration_old = true;
  } break;


  default:
    error->all(FC_FILE_LINE_FUNC, "undefined velocity-verlet integration type.");
  }
}

void Velocity_verlet::step_part_I () {
  // NOTHING
}

void Velocity_verlet::step_part_II () {

  auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;

  const auto psize = pos.size();

  switch (type) {


  case 1: {
    for (unsigned int i=0; i<psize; i++) { 
      // r(t+dt) = r(t) + v(t)*dt + 1/2 * a(t) * dt^2
      pos [i] += vel [i] * dt  +  0.5 * acc[i] * dt * dt;

      // v(t+dt/2) = v(t) + (a(t) ) * dt / 2
      vel [i] += 0.5 * dt * (acc [i]);

    }
  }
  break;

  case 2:{


    for (unsigned int i=0; i<psize; i++) { 
      // r(t+dt) = r(t) + v(t)*dt + 1/2 * a(t) * dt^2
      pos [i] += vel [i] * dt  +  0.5 * acc[i] * dt * dt;

    }
  } break;

  default:
    error->all(FC_FILE_LINE_FUNC, "undefined velocity-verlet integration type.");
  }

  
}


void Velocity_verlet::step_part_III () {

  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;
  auto &acc_old = atom_data -> owned.acceleration_old;

  switch (type) {



  case 1: {
    for (unsigned int i=0; i<vel.size(); i++) {
      // v(t+dt) = v(t+dt/2) + (a(t+dt) ) * dt / 2
      vel [i] += 0.5 * dt * (acc [i]);
    }
  }
  break;

  case 2:{
    for (unsigned int i=0; i<vel.size(); i++) {
      // v(t+dt) = v(t) + ( a(t+dt) + a(t) ) * dt / 2
      vel [i] += 0.5 * dt * (acc [i] + acc_old[i]);
    }
  } break;

  default:
    error->all(FC_FILE_LINE_FUNC, "undefined velocity-verlet integration type.");
  }
 
}

/*
FC_PYDEF_SETGET_PTR(Lj,atom_data,Atom_data);
FC_PYDEF_SETGET_PTR(Lj,domain,Domain);
FC_PYDEF_SETGET_PTR(Lj,neighborlist,Neighborlist);

FC_PYDEF_SETGET_STDVEC2D(Lj,epsilon,Real_t);  
FC_PYDEF_SETGET_STDVEC2D(Lj,sigma,Real_t);
FC_PYDEF_SETGET_STDVEC(Lj,epsilon_atom,Real_t);  
FC_PYDEF_SETGET_STDVEC(Lj,sigma_atom,Real_t);
FC_PYDEF_SETGET_STDVEC2D(Lj,cutoff_list,Real_t);

void export_py_Velocity_verlet () {

  using namespace boost::python;

  implicitly_convertible<std::shared_ptr<integrator::Velocity_verlet>,          
                         std::shared_ptr<Integrator> >(); 

  class_<integrator::Velocity_verlet>("Velocity_verlet",init<caviar::CAVIAR*>())
  ;

}
*/


} //integrator
} //objects

} // namespace caviar


