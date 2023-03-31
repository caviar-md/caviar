
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

#ifndef CAVIAR_OBJECTS_UNIQUE_MOLECULEGROUP_H
#define CAVIAR_OBJECTS_UNIQUE_MOLECULEGROUP_H

#include "caviar/objects/unique.h"
#include "caviar/objects/unique/molecule.h"

CAVIAR_NAMESPACE_OPEN


class Atom_data;
namespace unique {

/**
 * This class creates group of molecules.
 * groups create copies of the molecules.
 */
class Molecule_group : public Unique {
 public:
  Molecule_group (class CAVIAR *);
  Molecule_group (const Molecule_group &);
  Molecule_group ();
  ~Molecule_group ();

  bool read (caviar::interpreter::Parser *);
  void verify_settings ();

  Vector<double> pos_tot () const;
  Vector<double> vel_tot () const; 

  void add_molecule(const unique::Molecule &);
  void add_molecule(const unique::Molecule &,
                    caviar::Vector<double> p=caviar::Vector<double>{0,0,0},
                    caviar::Vector<double> v=caviar::Vector<double>{0,0,0});

  std::vector<unique::Molecule> molecules;

  bool part_of_a_molecule_group;    
  Molecule_group * upper_level_molecule_group;

  Vector<double> position, velocity;
 
};

} //unique


CAVIAR_NAMESPACE_CLOSE

#endif
