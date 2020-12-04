
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
FC_OBJECT_CREATOR_DEFAULT_FUNCTION(postprocess) {

  FC_GET_OBJECT_TYPE_AND_NAME

  objects::Postprocess * p_sh = nullptr; 

#include "caviar/objects/postprocess/macro/all.h"

  FC_ADD_OBJECT_TO_CONTAINER(postprocess)
}
#endif
