
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

#include "caviar/interpreter/object_handler/gdst.h"

CAVIAR_NAMESPACE_OPEN
namespace interpreter {
namespace object_handler {

int gdst (const std::string & t) { //get_dictionary_second_type 
  for (int i = 0; i <  static_cast<int>(object_handler::dictionary_second_type.size()); ++i) {
    if (dictionary_second_type[i]==t) return i;
  }
  return 0;
}

} //object_handler
} //interpreter

CAVIAR_NAMESPACE_CLOSE

