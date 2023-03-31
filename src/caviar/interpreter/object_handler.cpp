
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

#include "caviar/interpreter/object_handler.h"
#include "caviar/interpreter/object_handler/all.h"
#include "caviar/interpreter/object_handler/commands_map.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/all.h" // this is needed b.c. we want to call read()

CAVIAR_NAMESPACE_OPEN
namespace interpreter
{
  Object_handler::Object_handler(CAVIAR *fptr) : Pointers{fptr} {}

  Object_handler::~Object_handler()
  {
  }

  bool Object_handler::read_object(Parser *parser, const std::string object_name)
  {

    bool in_file = true;
    object_container = fptr->object_container;

    std::map<std::string, caviar::interpreter::object_handler::Dictionary>::iterator it;
    it = object_container->dictionary.find(object_name);

    if (it == object_container->dictionary.end())
      error->all(FC_FILE_LINE_FUNC_PARSE, static_cast<std::string>("Invalid object name: ") + object_name);

#define FC_LOCAL_CHECK_AND_RUN_READ_FUNCTION(OBJECT_TYPE)                    \
  if (it->second.type == object_handler::gdst(#OBJECT_TYPE))                 \
  {                                                                          \
    in_file = object_container->OBJECT_TYPE[it->second.index]->read(parser); \
    return in_file;                                                          \
  }

#define FC_GENERAL_CLASSNAME_MACRO(VAR1, VAR2, VAR3) \
  FC_LOCAL_CHECK_AND_RUN_READ_FUNCTION(VAR2)

#define FC_GENERAL_CLASSNAME_MACRO_ACTIVATED

#include "caviar/objects/macro/declaration/all.h"

#undef FC_GENERAL_CLASSNAME_MACRO_ACTIVATED
#undef FC_GENERAL_CLASSNAME_MACRO

#undef FC_LOCAL_CHECK_AND_RUN_READ_FUNCTION

    // basic types reads

    if (it->second.type == object_handler::gdst("boolean_variable"))
    {
      auto t = parser->get_raw_token();
      if (t.kind == caviar::interpreter::Kind::assign)
      {
        auto i = parser->get_bool();
        object_container->boolean_variable[it->second.index] = i;
        return true;
      }
      else
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "invalid operation");
      }
    }

    if (it->second.type == object_handler::gdst("int_variable"))
    {
      auto t = parser->get_raw_token();
      if (t.kind == caviar::interpreter::Kind::assign)
      {
        auto i = parser->get_int();
        object_container->int_variable[it->second.index] = i;
        return true;
      }
      else
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "invalid operation");
      }
    }

    if (it->second.type == object_handler::gdst("real_variable"))
    {
      auto t = parser->get_raw_token();
      if (t.kind == caviar::interpreter::Kind::assign)
      {
        auto i = parser->get_real();
        object_container->real_variable[it->second.index] = i;
        return true;
      }
      else
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "invalid operation");
      }
    }

    if (it->second.type == object_handler::gdst("real_3d_vector"))
    {
      auto t = parser->get_raw_token();
      if (t.kind == caviar::interpreter::Kind::assign)
      {
        auto ix = parser->get_real();
        auto iy = parser->get_real();
        auto iz = parser->get_real();
        object_container->real_3d_vector[it->second.index].x = ix;
        object_container->real_3d_vector[it->second.index].y = iy;
        object_container->real_3d_vector[it->second.index].z = iz;
        return true;
      }
      else
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "invalid operation");
      }
    }

    if (it->second.type == object_handler::gdst("int_3d_vector"))
    {
      auto t = parser->get_raw_token();
      if (t.kind == caviar::interpreter::Kind::assign)
      {
        auto ix = parser->get_int();
        auto iy = parser->get_int();
        auto iz = parser->get_int();
        object_container->int_3d_vector[it->second.index].x = ix;
        object_container->int_3d_vector[it->second.index].y = iy;
        object_container->int_3d_vector[it->second.index].z = iz;
        return true;
      }
      else
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "invalid operation");
      }
    }

    if (it->second.type == object_handler::gdst("string_variable"))
    {
      auto t = parser->get_raw_token();
      if (t.kind == caviar::interpreter::Kind::assign)
      {
        auto i = parser->get_string();
        object_container->string_variable[it->second.index] = i;
        return true;
      }
      else
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "invalid operation");
      }
    }

    // -----------

    return in_file; // WARNING
  }
} // interpreter
CAVIAR_NAMESPACE_CLOSE
