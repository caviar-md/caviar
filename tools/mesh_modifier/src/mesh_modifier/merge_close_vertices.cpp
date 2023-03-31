
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

#include "mesh_modifier.h"

#include <iostream>

namespace mesh_modifier
{

  void Mesh_modifier::merge_close_vertices(const double tolerance)
  {
    unsigned uci = unv_container.size() - 1;
    if (uci < 0)
    {
      std::cout << "error: there's no unv_container\n";
      return;
    }
  }
}
