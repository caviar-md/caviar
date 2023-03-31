
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



#ifdef FC_OBJECT_CREATOR_FUNCTION_DEFINITON
FC_OBJECT_CREATOR_DEFAULT_FUNCTION(force_field) {

  FC_GET_OBJECT_TYPE_AND_NAME

  Force_field * p_sh = nullptr; 

#include "caviar/objects/force_field/macro/all.h"

  FC_ADD_OBJECT_TO_CONTAINER(force_field)
}
#endif
