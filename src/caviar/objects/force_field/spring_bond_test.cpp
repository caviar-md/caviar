
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

#include "caviar/objects/force_field/spring_bond_test.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

namespace caviar {
namespace objects {
namespace force_field {

Spring_bond_test::Spring_bond_test (CAVIAR *fptr) : Force_field{fptr}
 {
  FC_OBJECT_INITIALIZE_INFO
}

bool Spring_bond_test::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"elastic_coef")) {
      GET_A_STDVECTOR_REAL_ELEMENT(elastic_coef)
      if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Elastic coef. have to be non-negative.");      
    }  else if (string_cmp(t,"dissip_coef")) {
      GET_A_STDVECTOR_REAL_ELEMENT(dissip_coef)
      if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Dissipation coef. have to be non-negative.");            
    } else if (string_cmp(t,"set_neighborlist") || string_cmp(t,"neighborlist")) {
      FIND_OBJECT_BY_NAME(neighborlist,it)
      neighborlist = object_container->neighborlist[it->second.index];
    } else if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"set_domain") || string_cmp(t,"domain")) {
      FIND_OBJECT_BY_NAME(domain,it)
      domain = object_container->domain[it->second.index];
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  
  return in_file;
}

void Spring_bond_test::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)
  FC_NULLPTR_CHECK(domain)
}

void Spring_bond_test::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS


  auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;
  auto &type = atom_data -> owned.type;
  auto &mass_inv = atom_data -> owned.mass_inv;

  auto &atomic_bond_index_vector = atom_data -> owned.atomic_bond_index_vector;
  auto &atomic_bond_vector = atom_data -> owned.atomic_bond_vector;


  for (unsigned int i=0; i<atomic_bond_index_vector.size(); i++) { 

    auto Nc = atomic_bond_index_vector[i].size();
    if (Nc==0) continue;
    for (unsigned int j=0; j<atomic_bond_vector[i].size(); j++) { 
        int k1 = atomic_bond_vector[i][j].index_1, k2 = atomic_bond_vector[i][j].index_2;
        int btype = atomic_bond_vector[i][j].type;
        double d = atomic_bond_vector[i][j].length;

#if defined(CAVIAR_WITH_MPI)
        const auto dr = pos[k2] - pos[k1]; 

#else    
        const auto dr = domain-> periodic_distance(pos[k2] - pos[k1]); 
#endif
        const auto dv = vel[k2] - vel[k1];
 
        const auto dr_sq = dr*dr;
        const auto dr_norm = std::sqrt(dr_sq);
        const auto dr_vec = dr / dr_norm;
        const auto force = -elastic_coef[btype]*(dr_norm - d)*dr_vec -(dissip_coef[btype] * dv);
//std::cout << force << std::endl;
        atom_data -> owned.acceleration [k1] -= force * mass_inv[type[k1]];
        atom_data -> owned.acceleration [k2] += force * mass_inv[type[k2]];        
    }

  }


}

} //force_field
} //objects
} // namespace caviar

