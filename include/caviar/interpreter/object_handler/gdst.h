
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

#ifndef CAVIAR_INTERPRETER_OBJECTHANDLER_GDST_H
#define CAVIAR_INTERPRETER_OBJECTHANDLER_GDST_H

#include "caviar/utility/caviar_config.h"
#include <string>
#include <vector>
//#include <map>

namespace caviar {
namespace interpreter {
namespace object_handler {


/**
 * This is a object to search and get a object index by means of its name.
 * 
 * 
 */
int gdst(const std::string &);//get_dictionary_second_type
//const static std::map<std::string,int> dictionary_second_type;

const std::vector<std::string> dictionary_second_type = {
  "wrong_type",

#define FC_GENERAL_CLASSNAME_MACRO(VAR1,VAR2,VAR3) \
  #VAR2 , 

#define FC_GENERAL_CLASSNAME_MACRO_ACTIVATED

#include "caviar/objects/macro/declaration/all.h"

#undef FC_GENERAL_CLASSNAME_MACRO_ACTIVATED
#undef FC_GENERAL_CLASSNAME_MACRO

#define FC_BASIC_TYPES_MACRO(VAR1)\
  #VAR1 ,

#define FC_BASIC_TYPES_MACRO_ACTIVATED
#include "caviar/interpreter/object_handler/all_basic_types_macro.h"
#undef FC_BASIC_TYPES_MACRO_ACTIVATED
#undef FC_BASIC_TYPES_MACRO
};

}
} //interpreter
} // namespace caviar

#endif
 
