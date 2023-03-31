
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

#include "caviar/objects/force_field/electrostatic_spherical_boundary.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"

#include <cmath>


CAVIAR_NAMESPACE_OPEN

namespace force_field {

double Electrostatic_spherical_boundary::potential (const Vector<double> &r) {

  initialize();

  double potential_sum = 0 ;  

  //const auto &pos = atom_data -> owned.position;  

  // particle-particle interaction part. If uncommented, it affects energy function.
  /* 
  for (unsigned int j=0;j<pos.size();++j) {

    const auto type_j = atom_data -> owned.type [j] ;
    const auto charge_j = atom_data -> owned.charge [ type_j ];      
    const auto dr = r - pos[j]; 
    const auto dr_sq = dr*dr;
    if (dr_sq == 0.0) continue;
    const auto dr_norm = std::sqrt(dr_sq);      
    potential_sum += charge_j  / dr_norm;
  }
  */

  // image-particle interaction
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:potential_sum)
#endif         
  for (unsigned int j=0;j<image.position.size();++j) {
    const auto charge_j = image.charge [ j ];      
    const auto dr = r - image.position[j]; 
    const auto dr_sq = dr*dr;
    const auto dr_norm = std::sqrt(dr_sq);      
    potential_sum += charge_j  / dr_norm;
  }

  return potential_sum * k_electrostatic;
}

double Electrostatic_spherical_boundary::potential (const int i) {

  initialize();

  double potential_sum = 0 ;  

  const auto &pos = atom_data -> owned.position;  

  // particle-particle interaction part. If uncommented, it affects energy function.
  /*
  for (unsigned int j=0;j<pos.size();++j) {
    if (i==static_cast<int>(j)) continue;
    const auto type_j = atom_data -> owned.type [j] ;
    const auto charge_j = atom_data -> owned.charge [ type_j ];      
    const auto dr = pos[i] - pos[j]; 
    const auto dr_sq = dr*dr;

    const auto dr_norm = std::sqrt(dr_sq);      
    potential_sum += charge_j  / dr_norm;
  }
  */

  // image-particle interaction
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:potential_sum)
#endif         
  for (unsigned int j=0;j<image.position.size();++j) {
    const auto charge_j = image.charge [ j ];      
    const auto dr = pos[i] - image.position[j]; 
    const auto dr_sq = dr*dr;
    const auto dr_norm = std::sqrt(dr_sq);      
    potential_sum += charge_j  / dr_norm;
  }

  return potential_sum * k_electrostatic;
}

} //force_field

CAVIAR_NAMESPACE_CLOSE

