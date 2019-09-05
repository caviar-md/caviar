
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

#ifndef CAVIAR_OBJECTS_ATOMDATA_UTILITY_ANGLE_H
#define CAVIAR_OBJECTS_ATOMDATA_UTILITY_ANGLE_H

namespace caviar {
namespace objects {
namespace atom_data {
// Angle contain data for rigid atomic angles which may be used in
// constraint algorithms or soft atomic angles in spring_angle force_fields
struct Angle {
  int index_1, index_2, index_3;  // atom index. 'index_2' is for the middle atom.
  int type; // used in soft atomic angles in force_fields
  double value; // angle value stored in radians. // TODO this can be stored by type in Atom_data
};
}
}
} // namespace caviar
#endif
