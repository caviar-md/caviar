
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

#include "caviar/objects/force_field/geometry.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/shape.h"
#include "caviar/objects/atom_data.h"
#include <string>
#include <cmath>
#include <fstream>

namespace caviar {
namespace objects {
namespace force_field {

Geometry::Geometry (CAVIAR *fptr) : Force_field {fptr}, 
shape_size_warning{false}
{ 
  FC_OBJECT_INITIALIZE_INFO
  young_modulus = 100.0;
  dissip_coef = 0.0;
}

Geometry::~Geometry () { 

}

bool Geometry::read (class caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"cutoff")) {
      error->all (FC_FILE_LINE_FUNC_PARSE, "cutoff is not implemented in force_field geometry");  
    } else if (string_cmp(t,"add_shape")) {
      FIND_OBJECT_BY_NAME(shape,it)
      shape.push_back (object_container->shape[it->second.index]);
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
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  
  return in_file;
}


void Geometry::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)

  if ((shape.size()==0) && (!shape_size_warning)) {
    output->warning("Geometry::calculate_acceleration: shape.size()==0");
    shape_size_warning = true;  
  }
}


void Geometry::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS


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
    //const auto r = radius [ type_i ];
     const auto r = a_radius [ type_i ];
    //Vector <Real_t> contact_vector {0,0,0};
    for (unsigned int j=0; j<shape.size();++j) {
      Vector <Real_t> contact_vector {0,0,0};        
      if (shape[j] -> in_contact(pos[i], r, contact_vector)) {
        acc[i] -= mass_inv_i * contact_vector*(young_modulus
                + (vel[i]*contact_vector)*dissip_coef/(contact_vector*contact_vector));
      }
    }
  }
}

} //force_field
} //objects
} // namespace caviar

