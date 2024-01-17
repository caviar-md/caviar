
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

#ifndef CAVIAR_INTERPRETER_OBJECTHANDLER_PREPROCESSORSNEW_H
#define CAVIAR_INTERPRETER_OBJECTHANDLER_PREPROCESSORSNEW_H

// This file have all it needs to make the code smaller, more generic and
// developer friendly. It contains macro functions that are used at the
// 'read (caviar::interpreter::Parser*)' function of the objects.

#define FC_IF_RAW_TOKEN_EOF_EOL                      \
  auto token = parser->get_raw_token();                  \
  std::string t = token.string_value;                   \
  if (token.kind == caviar::interpreter::Kind::eof)      \
  {                                                  \
    break;                                           \
  }                                                  \
  else if (token.kind == caviar::interpreter::Kind::eol) \
  {                                                  \
    break;                                           \
  }

// for now, we give the option of two different assignments :
// 'radius = 1 'radius 1'
// by this macro, we can change it all in the future
#define FC_HANDLE_ASSIGN                           \
  token = parser->get_raw_token();                     \
  if (token.kind == caviar::interpreter::Kind::assign) \
  {                                                \
  }                                                \
  else                                             \
    parser->keep_current_token();

#define FC_IF_GET_REAL(VARIABLE)   \
  if (string_cmp(t, #VARIABLE))   \
  {                                \
    FC_HANDLE_ASSIGN               \
    VARIABLE = parser->get_real(); \
    continue;                      \
  }

#define FC_IF_GET_INT(VARIABLE)   \
  if (string_cmp(t, #VARIABLE))  \
  {                               \
    FC_HANDLE_ASSIGN              \
    VARIABLE = parser->get_int(); \
    continue;                     \
  }

#define FC_IF_GET_BOOL(VARIABLE)   \
  if (string_cmp(t, #VARIABLE))   \
  {                                \
    FC_HANDLE_ASSIGN               \
    VARIABLE = parser->get_bool(); \
    continue;                      \
  }

#define FC_IF_GET_POSITIVE_REAL(VARIABLE)                                                                     \
  if (string_cmp(t, #VARIABLE))                                                                              \
  {                                                                                                           \
    FC_HANDLE_ASSIGN                                                                                          \
    VARIABLE = parser->get_real();                                                                            \
    if (VARIABLE <= 0)                                                                                        \
      error->all(FC_FILE_LINE_FUNC_PARSE, "expected a positive real '" + t + "' for object '" + OBJECT "'"); \
    continue;                                                                                                 \
  }

#define FC_IF_GET_POSITIVE_INT(VARIABLE)                                                                     \
  if (string_cmp(t, #VARIABLE))                                                                             \
  {                                                                                                          \
    FC_HANDLE_ASSIGN                                                                                         \
    VARIABLE = parser->get_int();                                                                            \
    if (VARIABLE <= 0)                                                                                       \
      error->all(FC_FILE_LINE_FUNC_PARSE, "expected a positive int '" + t + "' for object '" + OBJECT "'"); \
    continue;                                                                                                \
  }

#define FC_IF_GET_REAL3D(VARIABLE)   \
  if (string_cmp(t, #VARIABLE))     \
  {                                  \
    FC_HANDLE_ASSIGN                 \
    VARIABLE.x = parser->get_real(); \
    VARIABLE.y = parser->get_real(); \
    VARIABLE.z = parser->get_real(); \
    continue;                        \
  }

#define FC_IF_GET_INT3D(VARIABLE)   \
  if (string_cmp(t, #VARIABLE))    \
  {                                 \
    FC_HANDLE_ASSIGN                \
    VARIABLE.x = parser->get_int(); \
    VARIABLE.y = parser->get_int(); \
    VARIABLE.z = parser->get_int(); \
    continue;                       \
  }

#define FC_ERROR_PARAMETER(OBJECT) \
  error->all(FC_FILE_LINE_FUNC_PARSE, "unknown parameter '" + t + "' for object '" + OBJECT "'");

#endif
