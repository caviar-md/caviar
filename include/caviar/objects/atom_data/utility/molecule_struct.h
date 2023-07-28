
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

#ifndef CAVIAR_OBJECTS_ATOMDATA_UTILITY_MOLECULESTRUCT_H
#define CAVIAR_OBJECTS_ATOMDATA_UTILITY_MOLECULESTRUCT_H

#include "caviar/utility/objects_common_headers.h"
#include "caviar/objects/atom_data/utility/bond.h"
#include "caviar/objects/atom_data/utility/angle.h"
#include "caviar/objects/atom_data/utility/proper_dihedral.h"

CAVIAR_NAMESPACE_OPEN

namespace atom_data
{
 /**
   * It contains all the physical data of atoms and molecules
   * but since it is used only once here, its name won't be of any use. Also
   * a good name would be 'atom_data' which is used before!.
   */
  struct Molecule_struct  
  {

    /**
     * The first std::vector, is the molecule index. the inner data contain bonds.
     */
    std::vector<std::vector<atom_data::Bond>> atomic_bond_vector;


    /**
     * The first index, meaning
     * the first std::vector, is the molecule index. the inner data contain angles.
     */
    std::vector<std::vector<atom_data::Angle>> atomic_angle_vector;

    /**
     * The first index, meaning
     * the first std::vector, is the molecule index. the inner data contain Proper_dihedral.
     */
    std::vector<std::vector<atom_data::Proper_dihedral>> atomic_properdihedral_vector;

  };
}
CAVIAR_NAMESPACE_CLOSE
#endif
