
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

#ifndef CAVIAR_OBJECTS_UNIQUE_H
#define CAVIAR_OBJECTS_UNIQUE_H

#include "caviar/utility/objects_common_headers.h"

namespace caviar {

namespace objects {

/**
 * This class is the base class for all the unique objects.
 * These type of objects does not need a base class, but in order to have them
 * organized, we create them as the child of this class.
 */
class Unique : public Pointers {
 public:
  Unique (class CAVIAR *);

  virtual ~Unique ( );

  virtual bool read (class caviar::interpreter::Parser *);

  FC_BASE_OBJECT_COMMON_TOOLS
};

} //objects

} // namespace caviar

#endif
