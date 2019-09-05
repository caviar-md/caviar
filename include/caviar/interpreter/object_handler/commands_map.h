
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

#ifndef CAVIAR_INTERPRETER_OBJECTHANDLER_COMMANDSMAP_H
#define CAVIAR_INTERPRETER_OBJECTHANDLER_COMMANDSMAP_H

#include "caviar/interpreter/object_handler.h"

namespace caviar {
namespace interpreter {
//using CommandFunc_object_handler = bool (Object_handler::*) (Parser *); // a pointer to boolean function of ...

/**
 * A map between command names and the related functions.
 * 
 * 
 */
const std::map<std::string, CommandFunc_object_handler> Object_handler::commands_map = {

};
} //interpreter
} // namespace caviar

#endif
