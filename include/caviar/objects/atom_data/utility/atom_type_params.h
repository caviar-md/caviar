
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

#ifndef CAVIAR_OBJECTS_ATOMDATA_UTILITY_ATOMTYPEPARAMS_H
#define CAVIAR_OBJECTS_ATOMDATA_UTILITY_ATOMTYPEPARAMS_H

#include "caviar/utility/objects_common_headers.h"
//#include "caviar/objects/atom_data/utility/bond.h"
//#include "caviar/objects/atom_data/utility/angle.h"
//#include "caviar/objects/atom_data/utility/proper_dihedral.h"

CAVIAR_NAMESPACE_OPEN

namespace atom_data
{
 /**
   * It contains all the physical data of atoms types.
   */
  struct Atom_type_params  
  {    

    /**
     * 'mass' of an atom defined by the type. The mass may be used in
     * center-of-mass calculations and other functions. Do not depercate it.
     */
    std::vector<Real_t> mass;

    /**
     * simply the inverse value of 'mass' of an atom defined by the type.
     * since mass inverse is used in acceleration calculations.
     *
     */
    std::vector<Real_t> mass_inv;

    /**
     * 'charge' of an atom defined by the type.
     */
    std::vector<Real_t> charge;

    /**
     * 'radius' of an atom defined by the type. The user and the developers are
     * free to use this variable (for now!).
     */
    std::vector<Real_t> radius;
   


  };
}
CAVIAR_NAMESPACE_CLOSE
#endif
