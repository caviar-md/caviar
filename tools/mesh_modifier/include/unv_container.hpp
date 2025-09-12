
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

#ifndef UNV_CONTAINER_H
#define UNV_CONTAINER_H

#include <vector>
#include "universal_dataset_number_all.h"
#include "vector.h"

namespace mesh_modifier
{

  class Unv_container
  {
  public:
    Unv_container();
    ~Unv_container();

    std::vector<Universal_dataset_number_2411> udn_2411;
    std::vector<Universal_dataset_number_2412> udn_2412;
    std::vector<Universal_dataset_number_2467> udn_2467;
    std::vector<Universal_dataset_number_unsupported> udn_unsupported;
  };

}

#endif
