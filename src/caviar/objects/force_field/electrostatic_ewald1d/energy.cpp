
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

#include "caviar/objects/force_field/electrostatic_ewald1d.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"

#include <cmath>

namespace caviar {
namespace objects {
namespace force_field {

double Electrostatic_ewald1d::energy () {
// /* // XXX scheme using potential formula.
  const auto &pos = atom_data -> owned.position;    
  double energy_r = 0 ;
  for (unsigned int j=0;j<pos.size();++j) {
    const auto type_j = atom_data -> owned.type [j] ;  
    const auto charge_j = atom_data -> owned.charge [ type_j ];
//    energy_r += charge_j * potential(j); // 
    energy_r += charge_j * potential(pos [j]); //
  }
  return 0.5 * energy_r ;
// */

}

} //force_field
} //objects
} // namespace caviar

