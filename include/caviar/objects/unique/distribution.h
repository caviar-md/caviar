
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

#ifndef CAVIAR_OBJECTS_UNIQUE_DISTRIBUTION_H
#define CAVIAR_OBJECTS_UNIQUE_DISTRIBUTION_H

#include "caviar/objects/unique.h"

namespace caviar {


class Shape;
class Atom_data;
namespace unique {
class Molecule;
class Molecule_group;
class Atom;
class Atom_group;
class Grid_1D;
class Random_1D;

/**
 * This class creates initial arrangement of the particles using user's inputs.
 * It can put the atoms and molecules inside shapes.
 */
class Distribution : public Unique {
 public:
  Distribution (class CAVIAR *);
  ~Distribution () ;
  bool read (caviar::interpreter::Parser *);    
  void verify_settings ();
  bool distribute_grid_3D();
  bool distribute_random_3D(const int num, const double r);

  bool check_radius;
  
  class Atom_data *atom_data;
  class Shape *boundary_shape;
  class Atom *atom;
  class Atom_group *atom_group;
  class Molecule *molecule;
  class Molecule_group *molecule_group;

  class Grid_1D *grid_1d_x, *grid_1d_y, *grid_1d_z;
  class Random_1D *random_1d_x, *random_1d_y, *random_1d_z;
  
  std::vector<double> radius_vector;

};

} //unique


} // namespace caviar

#endif
