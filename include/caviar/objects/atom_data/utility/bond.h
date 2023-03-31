
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

#ifndef CAVIAR_OBJECTS_ATOMDATA_UTILITY_BOND_H
#define CAVIAR_OBJECTS_ATOMDATA_UTILITY_BOND_H

CAVIAR_NAMESPACE_OPEN

namespace atom_data {
// Bonds contain data for rigid atomic bonds which may be used in Shake like
// constraint algorithms or soft atomic bonds in spring_bond force_fields
struct Bond {
  int index_1, index_2;  // atom index
  int type; // used in soft atomic bonds in force_fields
  double length; // bond length // TODO this can be stored by type in Atom_data
};
}
CAVIAR_NAMESPACE_CLOSE
#endif
