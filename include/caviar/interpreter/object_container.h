
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

#ifndef CAVIAR_INTERPRETER_OBJECTCONTAINER_H
#define CAVIAR_INTERPRETER_OBJECTCONTAINER_H

#include "caviar/utility/pointers.h"
#include "caviar/utility/vector.h"
#include "caviar/utility/vector2D.h"
#include "caviar/interpreter/object_handler/dictionary.h"

#include <vector>
#include <string>
#include <map>
#include <unordered_set>
#include <memory>

#define FC_COMPLETE_FORWARD_DECLERATION
#include "caviar/objects/macro/declaration/all.h"
#undef FC_COMPLETE_FORWARD_DECLERATION

namespace caviar {
namespace interpreter {
class Parser;

/**
 * This class has the containers of all the objects.
 * 
 * 
 */
class Object_container : public Pointers  {
public:
  Object_container (class CAVIAR *);
  ~Object_container ();
  bool read (caviar::interpreter::Parser *);
  void report();

  std::map<std::string,caviar::interpreter::object_handler::Dictionary> dictionary;

  std::unordered_set<std::string> all_names;

#define FC_GENERAL_CLASSNAME_MACRO(VAR1,VAR2,VAR3) \
  std::vector< VAR3 *> VAR2; // 1

#define FC_GENERAL_CLASSNAME_MACRO_ACTIVATED

#include "caviar/objects/macro/declaration/all.h"

#undef FC_GENERAL_CLASSNAME_MACRO_ACTIVATED
#undef FC_GENERAL_CLASSNAME_MACRO

  // basic types
  std::vector<int> int_variable; // -1
  std::vector<double> real_variable; // -2
  std::vector<Vector2D<int>> int_2d_vector; // -3
  std::vector<Vector2D<double>> real_2d_vector; // -4
  std::vector<Vector<int>> int_3d_vector; // -5
  std::vector<Vector<double>> real_3d_vector; // -6
  std::vector<std::string> string_variable; // -7
  std::vector<bool> boolean_variable; // -8
  std::vector<std::shared_ptr<std::ofstream>>  ofs_objects;

public:

} ;
} //interpreter
} // namespace caviar

#endif
 
