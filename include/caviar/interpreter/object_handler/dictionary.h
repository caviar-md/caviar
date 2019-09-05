
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

#ifndef CAVIAR_INTERPRETER_OBJECTHANDLER_DICTIONARY_H
#define CAVIAR_INTERPRETER_OBJECTHANDLER_DICTIONARY_H

#include "caviar/utility/caviar_config.h"

namespace caviar {
namespace interpreter {
namespace object_handler {

/**
 * This class is a simple structure relating object types and index.
 * 
 * 
 */
class Dictionary {
  public:
  Dictionary () { };
  Dictionary (int t, int i) : type(t), index(i){ };
  ~Dictionary () {};
  int type; 
  int index; 

};

} // object_handler
} //interpreter
} // namespace caviar

#endif
 
