
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

#include "caviar/objects/md_simulator/basic.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/utility/time_utility.h"
#include "caviar/interpreter/communicator.h"
#ifdef CAVIAR_WITH_OPENMP
#include <omp.h>
#endif

#include <ctime>

namespace caviar {
namespace objects {
namespace md_simulator {

Basic::Basic (CAVIAR *fptr) : Md_simulator{fptr}
{
  FC_OBJECT_INITIALIZE_INFO
}

Basic::~Basic () {}


} //md_simulator
} //objects
} // namespace caviar

