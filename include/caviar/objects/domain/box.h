
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

#ifndef CAVIAR_OBJECTS_DOMAIN_BOX_H
#define CAVIAR_OBJECTS_DOMAIN_BOX_H

#include "caviar/objects/domain.h"

CAVIAR_NAMESPACE_OPEN

namespace domain
{

  /**
   * This class creates a 3D box for a simulation domain
   */
  class Box : public Domain
  {
  public:
    Box(class CAVIAR *);
    ~Box();
  public:
  };

} // domain

CAVIAR_NAMESPACE_CLOSE

#endif
