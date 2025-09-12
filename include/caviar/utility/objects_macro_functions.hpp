
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

#ifndef CAVIAR_UTILITY_MACRO_FUNCTIONS_H
#define CAVIAR_UTILITY_MACRO_FUNCTIONS_H

#define FC_SET_OBJECT_TYPE                                 \
  if (object_full_type == "")                              \
  {                                                        \
    object_full_type = __func__;                           \
  }                                                        \
  else                                                     \
  {                                                        \
    object_full_type = object_full_type + "::" + __func__; \
  }                                                        \
  object_type = __func__;

#define FC_OBJECT_INITIALIZE \
  FC_SET_OBJECT_TYPE         \
  settings_verified = false;

#define FC_OBJECT_INITIALIZE_INFO \
  FC_OBJECT_INITIALIZE            \
  output->info_create(object_full_type);

#define FC_OBJECT_READ_INFO \
  //output->info_read(object_full_type);


#define FC_OBJECT_READ_INFO_STR \
{\
  std::string s = "(Call) "+ object_full_type + ".read(): " + object_name + " " + t;\
  output->info(s);\
}
  

#define FC_BASE_OBJECT_COMMON_TOOLS                                   \
public:                                                               \
  std::string object_base_class_name, object_class_name, object_name; \
  std::string object_full_type, object_type;                          \
  bool settings_verified;                                             \
  virtual void verify_settings();

#define FC_OBJECT_VERIFY_SETTINGS \
  if (!settings_verified)         \
  {                               \
    verify_settings();            \
    settings_verified = true;     \
  }

#define FC_ERR_NOT_IMPLEMENTED \
  error->all(FC_FILE_LINE_FUNC, static_cast<std::string>("This feature is not completely implemented yet."));

#define FC_ERR_NOT_IMPLEMENTED_VAR(VAR) \
  error->all(FC_FILE_LINE_FUNC, static_cast<std::string>("This feature '") + VAR + "' is not completely implemented yet.");

#define FC_ERR_UNDEFINED \
  error->all(FC_FILE_LINE_FUNC, static_cast<std::string>("Undefined variable or command."));

#define FC_ERR_UNDEFINED_VAR(VAR) \
  error->all(FC_FILE_LINE_FUNC, static_cast<std::string>("Undefined variable or command: '") + VAR + "'.");

#define FC_NULLPTR_CHECK(OBJECT) \
  if (OBJECT == nullptr)         \
    error->all(FC_FILE_LINE_FUNC, static_cast<std::string>("Encountered 'nullptr' for '") + #OBJECT + "'.");

#define FC_CHECK_OBJECT_CLASS_NAME(VARVECTOR, VARIT, VARNAME)                        \
  auto ocn_st = object_container->VARVECTOR[VARIT->second.index]->object_class_name; \
  if (!string_cmp(ocn_st, #VARNAME))                                                 \
    error->all(FC_FILE_LINE_FUNC,                                                    \
               static_cast<std::string>("expected a '") + #VARNAME + "'object but got a '" + ocn_st + "' object.");

#endif
