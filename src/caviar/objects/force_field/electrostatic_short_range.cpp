
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

#include "caviar/objects/force_field/electrostatic_short_range.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

namespace caviar {
namespace objects {
namespace force_field {

Electrostatic_short_range::Electrostatic_short_range (CAVIAR *fptr) : Force_field{fptr}, 
    k_electrostatic{1.0}, external_field{Vector<double>{0,0,0}},
    beta{1.0}, initialized{false}
 {
  FC_OBJECT_INITIALIZE_INFO
    cutoff = -1.0;
}

bool Electrostatic_short_range::read (caviar::interpreter::Parser *parser) {
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
    } else if (string_cmp(t,"beta")) {
      GET_OR_CHOOSE_A_REAL(beta,"","")    
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

void Electrostatic_short_range::initialize () {
  initialized = true;
  C = 1.0/(beta * std::pow(cutoff, beta+1.0));
  D = - (beta + 1.0) / (beta * cutoff);
  cutoff_sq = cutoff * cutoff;
}

void Electrostatic_short_range::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)
  //FC_NULLPTR_CHECK(domain)
  FC_NULLPTR_CHECK(neighborlist)
  if (cutoff < 0.0) error->all (FC_FILE_LINE_FUNC, "Force field cutoff have to non-negative.");      
}


void Electrostatic_short_range::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS
  if (!initialized) initialize();

  const auto &pos = atom_data -> owned.position;  


  const auto &nlist = neighborlist -> neighlist;
  for (unsigned int i=0; i<nlist.size (); ++i) {
    const auto type_i = atom_data -> owned.type [i] ;
    const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];
    const auto charge_i = atom_data -> owned.charge [ type_i ];

    for (auto j : nlist[i]) {
      bool is_ghost = j >= nlist.size();
      Vector<Real_t> pos_j;
      Real_t type_j, mass_inv_j, charge_j;
      if (is_ghost) {
        j -= nlist.size ();
        pos_j = atom_data->ghost.position [j];
        type_j = atom_data->ghost.type [j];

      } else {
        pos_j = atom_data->owned.position [j];
        type_j = atom_data->owned.type [j];

      }
      charge_j = atom_data -> owned.charge [ type_j ];
      mass_inv_j = atom_data->owned.mass_inv [ type_j ];

      const auto dr = pos_j - pos[i]; 

      const auto dr_sq = dr*dr;
      if (dr_sq >= cutoff_sq) continue;
      const auto dr_norm = std::sqrt(dr_sq);      
      const auto force = k_electrostatic * charge_i * charge_j * dr / (dr_sq*dr_norm);
      const auto force_shifted = force * (1.0 - std::pow(dr_norm/cutoff, beta+1));
      atom_data -> owned.acceleration [i] -= force_shifted * mass_inv_i;

      if (!is_ghost) 
        atom_data -> owned.acceleration [j] += force_shifted * mass_inv_j;

      }
    
      const auto force = external_field * charge_i;
      atom_data -> owned.acceleration [i] += force * mass_inv_i;    
    }  



}

} //force_field
} //objects
} // namespace caviar

