
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

#ifndef CAVIAR_OBJECTS_MDSIMULATOR_BASIC_H
#define CAVIAR_OBJECTS_MDSIMULATOR_BASIC_H

#include "caviar/objects/md_simulator.h"

namespace caviar {

namespace md_simulator {

/**
 * This class has the basic class implementation for md_simulator
 * 
 * 
 */
class Basic : public Md_simulator {
 public:
  Basic (class CAVIAR *);
   ~Basic ( );

};

} //md_simulator

} // namespace caviar

#endif
