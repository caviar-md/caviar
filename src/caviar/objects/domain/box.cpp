
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

#include "caviar/objects/domain/box.h"
//#include "caviar/interpreter/communicator.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace domain
{

  Box::Box(CAVIAR *fptr) : Domain{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
  } 

  Box::~Box() {}

} // domain

CAVIAR_NAMESPACE_CLOSE
