
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

#include "caviar/objects/force_field/geometry_slab.h"
#include "caviar/utility/interpreter_io_headers.h"
//#include "caviar/objects/shape.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/unique/time_function_3d.h"
#include <string>
#include <cmath>
#include <fstream>

namespace caviar {

namespace force_field {

Geometry_slab::Geometry_slab (CAVIAR *fptr) : Force_field {fptr}
{ 
  FC_OBJECT_INITIALIZE_INFO
  young_modulus = 100.0;
  dissip_coef = 0.0;
  slab_direction = 0;
  symmetric = 0;
}

Geometry_slab::~Geometry_slab () { 

}

bool Geometry_slab::read (class caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"cutoff")) {
      error->all (FC_FILE_LINE_FUNC_PARSE, "cutoff is not implemented in force_field geometry");  
    } else if (string_cmp(t,"slab_direction")) {
      GET_OR_CHOOSE_A_INT(slab_direction,"","")
    } else if (string_cmp(t,"slab_position")) {
      GET_OR_CHOOSE_A_REAL(slab_position,"","")
    } else if (string_cmp(t,"asymmetric")) {
      symmetric = -1;
    } else if (string_cmp(t,"symmetric")) {
      symmetric = 1;
    } else if (string_cmp(t,"young_modulus")) {
      GET_OR_CHOOSE_A_REAL(young_modulus,"","")
    } else if (string_cmp(t,"dissip_coef")) {
      GET_OR_CHOOSE_A_REAL(dissip_coef,"","")
    } else if (string_cmp(t,"radius")) {
      GET_A_STDVECTOR_REAL_ELEMENT(radius)
      if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "radius have to be non-negative.");      
    } else if (string_cmp(t,"set_neighborlist") || string_cmp(t,"neighborlist")) {
      FIND_OBJECT_BY_NAME(neighborlist,it)
      neighborlist = object_container->neighborlist[it->second.index];
    } else if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"set_position_offset")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,time_function_3d)
      unique::Time_function_3d *a = dynamic_cast<unique::Time_function_3d *>(object_container->unique[it->second.index]);
      position_offset = a;
    } else if (string_cmp(t,"set_velocity_offset")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,time_function_3d)
      unique::Time_function_3d *a = dynamic_cast<unique::Time_function_3d *>(object_container->unique[it->second.index]);
      position_offset = a;
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  
  return in_file;
}


void Geometry_slab::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)

  if (slab_direction == 0)
    error->all(FC_FILE_LINE_FUNC, "slab_direction == 0");

  if (symmetric == 0)
    error->all(FC_FILE_LINE_FUNC, "symmetric == 0; expected symmetric(1) or asymmetric(-1)");

  if (slab_direction < 0 && symmetric==1)
    output->warning("Geometry_slab::verify_settings:: slab_direction < 0 && symmetric");
}


void Geometry_slab::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS
  Vector<double> p_o {0,0,0};
  if (position_offset != nullptr) p_o = position_offset->current_value;
  Vector<double> v_o {0,0,0};
  if (velocity_offset != nullptr) v_o = velocity_offset->current_value;

  const auto &pos = atom_data -> owned.position;
  const auto &vel = atom_data -> owned.velocity;  
  auto &acc = atom_data -> owned.acceleration;
  auto a_radius = atom_data->owned.radius;

  Vector<double> contact_vector {0, 0, 0};
  Vector<double> abs_contact_vector {0, 0, 0};

  auto abs_slab_direction = (slab_direction < 0 ? -slab_direction : slab_direction);

  switch (abs_slab_direction) {
  case +1 :
    contact_vector.x = slab_direction;
    abs_contact_vector.x = abs_slab_direction;
    break;

  case +2 :
    contact_vector.y = slab_direction;
    abs_contact_vector.y = abs_slab_direction;
    break;

  case +3 :
    contact_vector.z = slab_direction;
    abs_contact_vector.z = abs_slab_direction;
    break;
  }

  if (symmetric==1) {

#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
    for (unsigned int i=0;i<pos.size();++i) {
      const auto type_i = atom_data -> owned.type [i] ;
      const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];   
      double dx = 0;
      switch (slab_direction) {
  
      case -1 :
      case +1 :
        dx = pos[i].x - slab_position - p_o.x;
        break;

      case -2 :
      case +2 :
        dx = pos[i].y - slab_position - p_o.y;
        break;

      case -3 :
      case +3 :
        dx = pos[i].z - slab_position - p_o.z;
        break;
      }

      double dist = (dx < 0 ? -dx : dx);

      double compression = a_radius[ type_i ] - dist;

      if (compression > 0) {
        double dx_norm = dx / dist;
        acc[i] += mass_inv_i * abs_contact_vector* dx_norm *(young_modulus*compression
                - (vel[i]*abs_contact_vector)*dissip_coef);
      }
    }




  } else { //---------------------asymmetric-----------------------

#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
    for (unsigned int i=0;i<pos.size();++i) {
      const auto type_i = atom_data -> owned.type [i] ;
      const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];   
      double dx = 0;
      switch (slab_direction) {
  
      case -1 :
      case +1 :
        dx = pos[i].x - slab_position - p_o.x;
        break;

      case -2 :
      case +2 :
        dx = pos[i].y - slab_position - p_o.y;
        break;

      case -3 :
      case +3 :
        dx = pos[i].z - slab_position - p_o.z;
        break;
      }
      
  
      double compression = a_radius[ type_i ] - dx*slab_direction;

      if (compression > 0) {

        acc[i] += mass_inv_i * contact_vector*(young_modulus*compression
                - ((vel[i]-v_o)*abs_contact_vector)*dissip_coef);
      }
    }

  }

}

} //force_field

} // namespace caviar

