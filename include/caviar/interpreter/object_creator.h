
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

#ifndef CAVIAR_INTERPRETER_OBJECTCREATOR_H
#define CAVIAR_INTERPRETER_OBJECTCREATOR_H

#include "caviar/utility/pointers.h"

#include <string>
#include <map>

CAVIAR_NAMESPACE_OPEN
namespace interpreter {
class Parser;
using CommandFunc_object_creator = bool (Object_creator::*) (Parser *); // a pointer to boolean function of ...


/**
 * This class handles the object creations.
 * 
 * 
 */
class Object_creator : public Pointers  {
public:
  Object_creator (class CAVIAR *);
  ~Object_creator ();
  
  const static std::map<std::string, CommandFunc_object_creator> commands_map;
  
  // objects creator function declerations.

#define FC_GENERAL_CLASSNAME_MACRO(VAR1,VAR2,VAR3) \
  bool VAR2 (Parser *);

#define FC_GENERAL_CLASSNAME_MACRO_ACTIVATED

#include "caviar/objects/macro/declaration/all.h"

#undef FC_GENERAL_CLASSNAME_MACRO_ACTIVATED
#undef FC_GENERAL_CLASSNAME_MACRO


  // basic types creator function declerations.

#define FC_BASIC_TYPES_MACRO(VAR1) \
  bool VAR1 (Parser *);

#define FC_BASIC_TYPES_MACRO_ACTIVATED

#include "caviar/interpreter/object_handler/all_basic_types_macro.h"

#undef FC_BASIC_TYPES_MACRO_ACTIVATED
#undef FC_BASIC_TYPES_MACRO
  
public:

} ;
} //interpreter
CAVIAR_NAMESPACE_CLOSE

#endif
 
