
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

#include "caviar/objects/atom_data.h"
#include "caviar/interpreter/communicator.h"
#include "caviar/interpreter/error.h"
#include "caviar/objects/domain.h"
#include "caviar/objects/unique/atom.h"
#include "caviar/objects/unique/molecule.h"

#include <algorithm>

CAVIAR_NAMESPACE_OPEN

bool Atom_data::empty_of_atoms(const Vector<Real_t> p, double radius)
{
  double radius_sq = radius * radius;
  for (unsigned int i = 0; i < atom_struct_owned.position.size(); ++i)
  {
    auto dp = atom_struct_owned.position[i] - p;
    auto dp_sq = dp * dp;
    if (dp_sq < radius_sq)
      return false;
  }
  return true;
}

bool Atom_data::empty_of_atoms(const Vector<Real_t>, int)
{

  error->all(FC_FILE_LINE_FUNC, "not implemented");
  /*
    double radius_sq = radius * radius;
    for (int i = 0; i < atom_struct_owned.position.size(); ++i) {
      auto dp = atom_struct_owned.position[i] - p;
      auto dp_sq = dp*dp;
      if (dp_sq < radius_sq) return false;
    }
  */
  return true;
}

bool Atom_data::empty_of_atoms(unique::Atom &a)
{
  auto rad_a = atom_type_params.radius[a.type];
  for (unsigned int i = 0; i < atom_struct_owned.position.size(); ++i)
  {
    auto dp = atom_struct_owned.position[i] - a.pos_tot();
    auto dp_sq = dp * dp;
    auto rad_sum = atom_type_params.radius[atom_struct_owned.type[i]] + rad_a;
    if (dp_sq < rad_sum * rad_sum)
      return false;
  }
  return true;
}

bool Atom_data::empty_of_atoms(unique::Molecule &m)
{
  for (auto &&a : m.atoms)
  {
    if (!empty_of_atoms(a))
      return false;
  }
  return true;
}

CAVIAR_NAMESPACE_CLOSE
