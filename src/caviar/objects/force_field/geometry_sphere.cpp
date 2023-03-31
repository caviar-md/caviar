
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

#include "caviar/objects/force_field/geometry_sphere.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/unique/time_function_3d.h"
#include <string>
#include <cmath>
#include <fstream>

CAVIAR_NAMESPACE_OPEN

namespace force_field {

Geometry_sphere::Geometry_sphere (CAVIAR *fptr) : Force_field {fptr}
{ 
  FC_OBJECT_INITIALIZE_INFO
  center = caviar::Vector<double> (0,0,0);  
  young_modulus = 100.0;
  dissip_coef = 0.0;
  inside = true;
  radius = 1.0;
}

Geometry_sphere::~Geometry_sphere () { 

}

bool Geometry_sphere::read (class caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"center")) {
      GET_OR_CHOOSE_A_REAL_3D_VECTOR(center,"","")
    } else if (string_cmp(t,"radius")) {
      GET_OR_CHOOSE_A_REAL(radius,"","")
      if (radius < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "radius have to be non-negative.");      
    } else if (string_cmp(t,"outside")) {
      inside = false;
    } else if (string_cmp(t,"inside")) {
      inside = true;
    } else if (string_cmp(t,"young_modulus")) {
      GET_OR_CHOOSE_A_REAL(young_modulus,"","")
    } else if (string_cmp(t,"dissip_coef")) {
      GET_OR_CHOOSE_A_REAL(dissip_coef,"","")
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


void Geometry_sphere::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)  
}


void Geometry_sphere::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS
  Vector<double> p_o {0,0,0};
  if (position_offset != nullptr) p_o = position_offset->current_value;
  Vector<double> v_o {0,0,0};
  if (velocity_offset != nullptr) v_o = velocity_offset->current_value;

  const auto &pos = atom_data -> owned.position;
  const auto &vel = atom_data -> owned.velocity;  
  auto &acc = atom_data -> owned.acceleration;
  auto a_radius = atom_data->owned.radius;

  
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
    for (unsigned int i=0;i<pos.size();++i) {
      const auto type_i = atom_data -> owned.type [i] ;
      const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];   
      
      auto cp = (pos[i]-center - p_o);
      auto cp_sq = cp*cp;
      auto cp_abs = std::sqrt(cp_sq);
      auto cp_norm = cp / cp_abs;
      auto vel_coef = 1.0; // TODO: CHECK THIS
      
      double compression;
      
      if (inside)
      {
          compression = cp_abs + a_radius [type_i] - radius;                      
          cp_norm *= -1.0;
          
      }
      else
      {
          compression = radius - (cp_abs - a_radius [type_i]) ; 
          //vel_coef = -1.0; // TODO: CHECK THIS for inside case
      }

      if (compression > 0) {

        acc[i] += mass_inv_i * cp_norm*(young_modulus*compression
                - ((vel[i]-v_o)*cp_norm)*dissip_coef*vel_coef);
        
      }
    }

}

} //force_field

CAVIAR_NAMESPACE_CLOSE

