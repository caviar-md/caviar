
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

#ifndef CAVIAR_OBJECTS_UNIQUE_ATOM_H
#define CAVIAR_OBJECTS_UNIQUE_ATOM_H

#include "caviar/objects/unique.h"
#include "caviar/utility/vector.h"
#include <vector>

namespace caviar {

namespace objects {
namespace unique {
class Molecule;
class Atom_group;


/**
 * This class creates atoms as a tool for initial position of the particles for
 *  atom_data class
 */
class Atom  : public Unique {
 public:
  Atom (class CAVIAR *);    
  Atom (const Atom &);
  Atom ();
  ~Atom () ;
  bool read (caviar::interpreter::Parser *);
  void verify_settings ();
  Vector<double> pos_tot () const;
  Vector<double> vel_tot () const; 

  void output_xyz (std::ofstream &);  
  void extract_all_e_pos_vel (std::vector<int>&, std::vector<Vector<double>>&, std::vector<Vector<double>>&);      

  bool part_of_a_molecule;    
  Molecule * upper_level_molecule;

  bool part_of_a_atom_group;    
  Atom_group * upper_level_atom_group;
  
  Vector<double> position, velocity;
  unsigned int type;

};

} //unique
} //objects

} // namespace caviar

#endif
