
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

#ifndef CAVIAR_OBJECTS_UNIQUE_MOLECULE_H
#define CAVIAR_OBJECTS_UNIQUE_MOLECULE_H

#include "caviar/objects/unique.h"
#include "caviar/objects/atom_data/utility/bond.h"
#include "caviar/objects/atom_data/utility/angle.h"
#include "caviar/objects/atom_data/utility/proper_dihedral.h"

CAVIAR_NAMESPACE_OPEN


namespace unique {
class Atom;
class Molecule_group;

/**
 * This class creates molecules as a tool for initial position of the particles for
 *  atom_data class
 */
class Molecule : public Unique {
 public:
  Molecule (class CAVIAR *);
  Molecule (const Molecule & a);
  Molecule ();
  ~Molecule ();

  bool read (caviar::interpreter::Parser *);
  void verify_settings ();

  Vector<double> pos_tot () const;
  Vector<double> vel_tot () const;  

  bool add_atom (caviar::interpreter::Parser *); 
  bool add_atom (const class Atom &);
  bool add_atom (const class Atom &, const Vector<double> &p, const Vector<double> &v);  

  // called by the molecule itself. The output file name is automaticly generated
  void output_xyz ();

  // could be called by a Molecule_group or Molecule_list
  void output_xyz (std::ofstream &);
  
  // could be called in a library-type call. get's the output file name.
  void output_xyz (const std::string &);

  // puts the atoms type, total position and velocity inside the vectors.
  void extract_all_e_pos_vel (std::vector<int>&, std::vector<Vector<double>>&,
      std::vector<Vector<double>>&);
  
  bool part_of_a_molecule_group;    
  Molecule_group * upper_level_molecule_group;

  Vector<double> position, velocity;
  std::vector<Atom> atoms;

  std::vector<atom_data::Bond> atomic_bond; 
  std::vector<int> atomic_bond_index; 

  std::vector<atom_data::Angle> atomic_angle; 
  std::vector<int> atomic_angle_index; 

  std::vector<atom_data::Proper_dihedral> atomic_properdihedral; 
  std::vector<int> atomic_properdihedral_index;

  
};

} //unique


CAVIAR_NAMESPACE_CLOSE

#endif
