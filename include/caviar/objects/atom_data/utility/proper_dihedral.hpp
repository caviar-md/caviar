
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

#ifndef CAVIAR_OBJECTS_ATOMDATA_UTILITY_PROPERDIHEDRAL_H
#define CAVIAR_OBJECTS_ATOMDATA_UTILITY_PROPERDIHEDRAL_H

CAVIAR_NAMESPACE_OPEN

namespace atom_data
{
  // Proper Dihedral contain data for rigid atomic proper dihedral which may be used in
  // constraint algorithms or soft atomic proper dihedral in harmonic_proper_dihedral force_fields
  struct Proper_dihedral
  {
    int id_1, id_2, id_3, id_4;
    int type; // used in soft atomic proper dihedral in force_fields
  };
}
CAVIAR_NAMESPACE_CLOSE
#endif
