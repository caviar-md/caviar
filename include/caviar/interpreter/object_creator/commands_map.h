
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

#ifndef CAVIAR_INTERPRETER_OBJECTCREATOR_COMMANDSMAP_H
#define CAVIAR_INTERPRETER_OBJECTCREATOR_COMMANDSMAP_H

#include "caviar/interpreter/object_creator.h"

CAVIAR_NAMESPACE_OPEN
namespace interpreter
{
    // using CommandFunc_object_creator = bool (Object_creator::*) (Parser *); // a pointer to boolean function of ...

    /**
     * A map between command names and the related functions.
     *
     *
     */
    const std::map<std::string, CommandFunc_object_creator> Object_creator::commands_map = {

#define FC_GENERAL_CLASSNAME_MACRO(VAR1, VAR2, VAR3) \
    {#VAR2, &Object_creator::VAR2},

#define FC_GENERAL_CLASSNAME_MACRO_ACTIVATED

#include "caviar/objects/macro/declaration/all.h"

#undef FC_GENERAL_CLASSNAME_MACRO_ACTIVATED
#undef FC_GENERAL_CLASSNAME_MACRO

        // basic types
        //{  "int_variable", &Object_creator::int_variable},
        {"int", &Object_creator::int_variable},
        //{  "real_variable", &Object_creator::real_variable},
        {"real", &Object_creator::real_variable},
        //{  "int_2d_vector", &Object_creator::int_2d_vector},
        //{  "int_2d", &Object_creator::int_2d_vector},
        {"int2d", &Object_creator::int_2d_vector},
        //{  "real_2d_vector", &Object_creator::real_2d_vector},
        //{  "real_2d", &Object_creator::real_2d_vector},
        {"real2d", &Object_creator::real_2d_vector},
        //{  "int_3d_vector", &Object_creator::int_3d_vector},
        //{  "int_3d", &Object_creator::int_3d_vector},
        {"int3d", &Object_creator::int_3d_vector},
        //{  "real_3d_vector", &Object_creator::real_3d_vector},
        //{  "real_3d", &Object_creator::real_3d_vector},
        {"real3d", &Object_creator::real_3d_vector},
        //{  "string_variable", &Object_creator::string_variable},
        {"string", &Object_creator::string_variable},
        //{  "boolean_variable", &Object_creator::boolean_variable},
        //{  "boolean", &Object_creator::boolean_variable},
        {"bool", &Object_creator::boolean_variable},
    };
} // interpreter
CAVIAR_NAMESPACE_CLOSE

#endif
