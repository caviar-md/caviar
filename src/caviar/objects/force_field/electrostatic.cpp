
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

#include "caviar/objects/force_field/electrostatic.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field {

Electrostatic::Electrostatic (CAVIAR *fptr) : Force_field{fptr}, k_electrostatic{1.0}, external_field{Vector<double>{0,0,0}}
 {
  FC_OBJECT_INITIALIZE_INFO
}

bool Electrostatic::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"cutoff")) {
      GET_OR_CHOOSE_A_REAL(cutoff,"","")
      if (cutoff < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "Force field cutoff have to non-negative.");      
    } else if (string_cmp(t,"external_field")) {
      GET_OR_CHOOSE_A_REAL_3D_VECTOR(external_field, "", "");
      //if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Epsilon have to be non-negative.");      
    }  else if (string_cmp(t,"k_electrostatic")) {
      GET_OR_CHOOSE_A_REAL(k_electrostatic,"","")    
      if (k_electrostatic < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "k_electrostatic has to be non-negative.");            
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

void Electrostatic::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)
  FC_NULLPTR_CHECK(domain)
}


void Electrostatic::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS


  const auto &pos = atom_data -> owned.position;  
  {
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
    for (unsigned int i=0;i<pos.size();++i) {
      const auto type_i = atom_data -> owned.type [i] ;

      const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];
      const auto charge_i = atom_data -> owned.charge [ type_i ];
      for (unsigned int j=i+1;j<pos.size();++j) {

        const auto type_j = atom_data -> owned.type [j] ;

        const auto mass_inv_j = atom_data -> owned.mass_inv [ type_j ];
        const auto charge_j = atom_data -> owned.charge [ type_j ];

#if defined(CAVIAR_WITH_MPI)
        const auto dr = pos[j] - pos[i]; 
#else     
        const auto dr = domain-> periodic_distance(pos[j] - pos[i]); 
#endif

        const auto dr_sq = dr*dr;
        const auto dr_norm = std::sqrt(dr_sq);      
        const auto force = k_electrostatic * charge_i * charge_j * dr / (dr_sq*dr_norm);
        atom_data -> owned.acceleration [i] -= force * mass_inv_i;
        
#ifdef CAVIAR_WITH_OPENMP        
#pragma omp atomic 
          atom_data -> owned.acceleration [j].x -= force.x * mass_inv_j;   
#pragma omp atomic
          atom_data -> owned.acceleration [j].y -= force.y * mass_inv_j;   
#pragma omp atomic 
          atom_data -> owned.acceleration [j].z -= force.z * mass_inv_j;  
#else
          atom_data -> owned.acceleration [j] -= force * mass_inv_j;   
#endif                        

      }
    
      const auto force = external_field * charge_i;
      atom_data -> owned.acceleration [i] += force * mass_inv_i;    
    }  

  }

}

} //force_field

CAVIAR_NAMESPACE_CLOSE

