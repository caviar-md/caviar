
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

#include "caviar/objects/force_field/electrostatic_ewald_r.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field {

double Electrostatic_ewald_r::energy () {
// /* // XXX working scheme using potential formula. order (neighborlist)
  const auto &pos = atom_data -> owned.position;    
  double energy_r = 0 ;
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:energy_r)
#endif    
  for (unsigned int j=0;j<pos.size();++j) {
    const auto type_j = atom_data -> owned.type [j] ;  
    const auto charge_j = atom_data -> owned.charge [ type_j ];
//    energy_r += charge_j * potential(j); // (neighlist) verlet_list and cell_list
    energy_r += charge_j * potential(pos [j]); // working (binlist) cell_list.
  }
  return 0.5 * energy_r ;
// */

/* // XXX working scheme with one sum on neighborlist

  double e = 0 ;  
  const auto &pos = atom_data -> owned.position;
  const auto &nlist = neighborlist -> neighlist;

  const unsigned pos_size = pos.size();

  for (unsigned int i=0; i<nlist.size (); ++i) {
    const auto &pos_i = atom_data -> owned.position [i];
    const auto type_i = atom_data -> owned.type [i];
    const auto charge_i = atom_data -> owned.charge [ type_i ]; 

    for (auto j : nlist[i]) {
      double coef = 2.0;
      bool is_ghost = j >= pos_size;

      Vector<Real_t> pos_j;
      Real_t type_j;
      if (is_ghost) {
        coef = 1.0;
        j -= pos_size;
        pos_j = atom_data->ghost.position [j];
        type_j = atom_data->ghost.type [j];
      } else {
        pos_j = atom_data->owned.position [j];
        type_j = atom_data->owned.type [j];
      }

      const auto charge_j = atom_data -> owned.charge [ type_j ]; 
      const auto r_ij = pos_i - pos_j;

      if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0) continue;
      
      const auto r_ij_norm = std::sqrt(r_ij*r_ij);
      const auto erfc_arg = alpha*r_ij_norm;

      e += coef * charge_i*charge_j * std::erfc(erfc_arg) / r_ij_norm;    
     
    }
  }

  return 0.5*e* k_electrostatic;
// */

/* //  XXX workin scheme without neighborlist

  double e = 0 ;  
  const auto &pos = atom_data -> owned.position;
  const auto &gpos = atom_data -> ghost.position;
  const unsigned pos_size = pos.size();

  for (unsigned int i=0; i<pos_size; ++i) {
    const auto &pos_i = atom_data -> owned.position [i];
    const auto type_i = atom_data -> owned.type [i];
    const auto charge_i = atom_data -> owned.charge [ type_i ]; 

    for (unsigned int j=0; j<pos_size; ++j) {
      const auto pos_j = atom_data->owned.position [j];
      const auto type_j = atom_data->owned.type [j];

      const auto charge_j = atom_data -> owned.charge [ type_j ]; 
      const auto r_ij = pos_i - pos_j;

      if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0) continue;

      const auto r_ij_norm = std::sqrt(r_ij*r_ij);
      const auto erfc_arg = alpha*r_ij_norm;

      e += charge_i*charge_j * std::erfc(erfc_arg) / r_ij_norm;    
    }

    for (unsigned int j=0; j<gpos.size(); ++j) {
      const auto pos_j = atom_data->ghost.position [j];
      const auto type_j = atom_data->ghost.type [j];

      const auto charge_j = atom_data -> owned.charge [ type_j ]; 
      const auto r_ij = pos_i - pos_j;

      //if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0) continue;

      const auto r_ij_norm = std::sqrt(r_ij*r_ij);
      const auto erfc_arg = alpha*r_ij_norm;

      e += charge_i*charge_j * std::erfc(erfc_arg) / r_ij_norm;    
    }

  }

  return 0.5*e* k_electrostatic;
*/

/* // XXX working scheme with neighborlist and two sums
  double e = 0 ;  
  const auto &pos = atom_data -> owned.position;
  const auto &gpos = atom_data -> ghost.position;
  const unsigned pos_size = pos.size();

  const auto &nlist = neighborlist -> neighlist;

  for (unsigned int i=0; i<pos_size; ++i) {
    const auto &pos_i = atom_data -> owned.position [i];
    const auto type_i = atom_data -> owned.type [i];
    const auto charge_i = atom_data -> owned.charge [ type_i ]; 

    double se1 = 0;

    for (auto j : nlist[i]) {
      bool is_ghost = j >= pos_size;
      if (is_ghost) continue;
      const auto pos_j = atom_data->owned.position [j];
      const auto type_j = atom_data->owned.type [j];

      const auto charge_j = atom_data -> owned.charge [ type_j ]; 
      const auto r_ij = pos_i - pos_j;

      if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0) continue;

      const auto r_ij_norm = std::sqrt(r_ij*r_ij);
      const auto erfc_arg = alpha*r_ij_norm;

      se1 += charge_i*charge_j * std::erfc(erfc_arg) / r_ij_norm;    
    }
    se1 *= 2.0;

    double se2 = 0;
    for (auto j : nlist[i]) {
      bool is_ghost = j >= pos_size;
      if (!is_ghost) continue;
      j -= pos_size;
      const auto pos_j = atom_data->ghost.position [j];
      const auto type_j = atom_data->ghost.type [j];

      const auto charge_j = atom_data -> owned.charge [ type_j ]; 
      const auto r_ij = pos_i - pos_j;

      //if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0) continue;

      const auto r_ij_norm = std::sqrt(r_ij*r_ij);
      const auto erfc_arg = alpha*r_ij_norm;

      se2 += charge_i*charge_j * std::erfc(erfc_arg) / r_ij_norm;    
    }
    e += se1 + se2;
  }

  return 0.5*e* k_electrostatic;
*/

}

} //force_field

CAVIAR_NAMESPACE_CLOSE

