
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

#ifndef CAVIAR_INTERPRETER_OBJECTHANDLER_H
#define CAVIAR_INTERPRETER_OBJECTHANDLER_H

#include "caviar/utility/pointers.h"
#include <string>
#include <map>

CAVIAR_NAMESPACE_OPEN
namespace interpreter {
class Parser;
using CommandFunc_object_handler = bool (Object_handler::*) (class caviar::interpreter::Parser *); // a pointer to boolean function of ...


/**
 * This class is a way to call all the created objects.
 * 
 * 
 */
class Object_handler : public Pointers  {
public:
  Object_handler (class CAVIAR *);
  ~Object_handler ();
  
  const static std::map<std::string, CommandFunc_object_handler> commands_map;
  
  bool read_object (Parser *, const std::string);

public:

} ;
} //interpreter
CAVIAR_NAMESPACE_CLOSE

#endif
 
