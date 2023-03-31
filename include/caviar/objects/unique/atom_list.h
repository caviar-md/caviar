
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

#ifndef CAVIAR_OBJECTS_UNIQUE_ATOMLIST_H
#define CAVIAR_OBJECTS_UNIQUE_ATOMLIST_H

#include "caviar/objects/unique.h"
#include "caviar/objects/unique/atom.h"

CAVIAR_NAMESPACE_OPEN


namespace unique {

/**
 * This class creates list of atoms.
 * list contains references of the atoms.
 */
class Atom_list  : public Unique {
 public:
  Atom_list (class CAVIAR *);    
  Atom_list (const Atom_list &);
  Atom_list ();
  ~Atom_list () ;
  bool read (caviar::interpreter::Parser *);
  void verify_settings ();
  void add_atom(const unique::Atom &);
  std::vector<unique::Atom *> atoms;

};

} //unique


CAVIAR_NAMESPACE_CLOSE

#endif
