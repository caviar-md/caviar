
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

#include "caviar/interpreter/object_creator.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/all_derived_classes.h"

CAVIAR_NAMESPACE_OPEN
namespace interpreter
{

#define FC_OBJECT_CREATOR_DEFAULT_FUNCTION(VAR1) \
  bool Object_creator::VAR1(Parser *parser)

#define FC_GET_OBJECT_TYPE_AND_NAME            \
  auto object_type = parser->get_identifier(); \
  auto object_name = parser->get_identifier(); \
  bool object_type_found = false;

#define FC_CHECK_AND_CREATE(VAR1, VAR2, VAR3)                   \
  if ((!object_type_found) && string_cmp_i(object_type, #VAR2)) \
  {                                                             \
    object_type_found = true;                                   \
    p_sh = new VAR3(fptr);                                      \
  }

#define FC_ADD_OBJECT_TO_CONTAINER(OBJECT_TYPE)                                                                              \
  if (!object_type_found)                                                                                                    \
    error->all(FC_FILE_LINE_FUNC_PARSE, static_cast<std::string>("Undefined '") + #OBJECT_TYPE + "::" + object_type + "'."); \
  p_sh->object_name = object_name;                                                                                           \
  p_sh->object_class_name = object_type;                                                                                     \
  p_sh->object_base_class_name = __func__;                                                                                   \
  object_container->all_names.insert(object_name);                                                                           \
  int index = object_container->OBJECT_TYPE.size();                                                                          \
  object_container->OBJECT_TYPE.emplace_back(p_sh);                                                                          \
  object_handler::Dictionary dict(object_handler::gdst(#OBJECT_TYPE), index);                                                \
  object_container->dictionary.insert(std::make_pair(object_name, dict));                                                    \
  return true;

#define FC_OBJECT_CREATOR_FUNCTION_DEFINITON
#define FC_CHECK_AND_CREATE_ACTIVATED

#include "caviar/objects/macro/creation/all.h"

#undef FC_CHECK_AND_CREATE_ACTIVATED
#undef FC_OBJECT_CREATOR_FUNCTION_DEFINITON
#undef FC_OBJECT_CREATOR_DEFAULT_FUNCTION
#undef FC_CHECK_AND_CREATE
#undef FC_GET_OBJECT_TYPE_AND_NAME
#undef FC_ADD_OBJECT_TO_CONTAINER

} // interpreter
CAVIAR_NAMESPACE_CLOSE
