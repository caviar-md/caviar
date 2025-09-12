
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

#ifndef CAVIAR_OBJECTS_UNIQUE_ATOMGROUP_H
#define CAVIAR_OBJECTS_UNIQUE_ATOMGROUP_H

#include "caviar/objects/unique.h"
#include "caviar/objects/unique/atom.h"

CAVIAR_NAMESPACE_OPEN

namespace unique
{

  /**
   * This class creates group of atoms.
   * groups create copies of the atoms.
   */
  class Atom_group : public Unique
  {
  public:
    Atom_group(class CAVIAR *);
    Atom_group(const Atom_group &);
    Atom_group();
    ~Atom_group();
    bool read(caviar::interpreter::Parser *);
    Vector<double> pos_tot() const;
    Vector<double> vel_tot() const;
    void add_atom(const unique::Atom &);
    void add_atom(const unique::Atom &,
                  caviar::Vector<double> p = caviar::Vector<double>{0, 0, 0},
                  caviar::Vector<double> v = caviar::Vector<double>{0, 0, 0});

    std::vector<unique::Atom> atoms;

    bool part_of_a_atom_group;
    Atom_group *upper_level_atom_group;

    Vector<double> position, velocity;
  };

} // unique

CAVIAR_NAMESPACE_CLOSE

#endif
