
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

#ifndef CAVIAR_OBJECTS_UNIQUE_MOLECULELIST_H
#define CAVIAR_OBJECTS_UNIQUE_MOLECULELIST_H

#include "caviar/objects/unique.h"
#include "caviar/objects/unique/molecule.h"

namespace caviar {


class Atom_data;
namespace unique {

/**
 * This class creates list of molecules.
 * list contains references of the molecules.
 */
class Molecule_list : public Unique {
 public:
  Molecule_list (class CAVIAR *);
  Molecule_list (const Molecule_list &);
  Molecule_list ();
  ~Molecule_list ();

  bool read (caviar::interpreter::Parser *);
  void verify_settings ();
  void add_molecule(const unique::Molecule &);
  std::vector<unique::Molecule*> molecules;

  
};

} //unique


} // namespace caviar

#endif
