
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

double Electrostatic_ewald_r::potential (const Vector<double> &r) {
// XXX Working scheme using binlist
  double potential_r = 0 ;  

  const auto &pos = atom_data -> owned.position;
  const auto &binlist = neighborlist -> binlist;
  const auto &nb = neighborlist -> neigh_bin;
  const auto nb_i = neighborlist -> neigh_bin_index (r);
  const int pos_size = pos.size();

  const auto pos_i = r;
  
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:potential_r)
#endif 
  for (unsigned nb_j = 0; nb_j < nb[nb_i].size(); ++nb_j) {
    const auto &nb_ij = nb[nb_i][nb_j];

    for (unsigned i = 0; i < binlist [nb_ij.x] [nb_ij.y] [nb_ij.z].size(); ++i) {

      int j = binlist[nb_ij.x] [nb_ij.y] [nb_ij.z][i];

      bool is_ghost = j >= pos_size;

      Vector<Real_t> pos_j;
      Real_t type_j;
      if (is_ghost) {

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

      potential_r += charge_j * std::erfc(erfc_arg) / r_ij_norm;    
    }
  }
  return potential_r * k_electrostatic;
}

double Electrostatic_ewald_r::potential (const int i) {
// XXX Working scheme (FOR ENERGY) using neighlist
  error->all("not implemented. It gives a correct answer for the energy but the formula is wrong for one point");

  double potential_r = 0 ;  
  const auto &pos = atom_data -> owned.position;
  const auto &nlist = neighborlist -> neighlist;

  const unsigned pos_size = pos.size();

  const auto &pos_i = atom_data -> owned.position [i];
  
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:potential_r)
#endif   
  for (unsigned int k = 0; k < nlist[i].size(); ++k) {
      auto j = nlist[i][k];
    
      double coef = 2.0; // ewald: 'coef=2' for owned in 'neighlist'. Not for binlist.
      bool is_ghost = j >= pos_size;

      Vector<Real_t> pos_j;
      Real_t type_j;
      if (is_ghost) {
        coef = 1.0; // ewald:'coef=1' for ghost in 'neighlist'. Not for binlist.
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

      potential_r += coef * charge_j * std::erfc(erfc_arg) / r_ij_norm;    
  }

  return potential_r * k_electrostatic;
}

} //force_field

CAVIAR_NAMESPACE_CLOSE

