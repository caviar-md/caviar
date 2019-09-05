
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

#ifndef CAVIAR_OBJECTS_MDSIMULATOR_POLYATOMIC_H
#define CAVIAR_OBJECTS_MDSIMULATOR_POLYATOMIC_H

#include "caviar/objects/md_simulator.h"

namespace caviar {
namespace objects {
namespace md_simulator {

/**
 * This class has an implementation for md_simulator according to
 * Berendsen et al, 1984 : 
 * J. Chem. Phys. 81, 3684 (1984); doi: 10.1063/1.448118
 * 
 */
class Polyatomic : public Md_simulator {
 public:
  Polyatomic (class CAVIAR *);
   ~Polyatomic ( );
  bool read (class caviar::interpreter::Parser *);
  bool run ();
  void verify_settings ();
  void step(int i);


};

} //md_simulator
} //objects
} // namespace caviar

#endif
