
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

#include "caviar/objects/constraint/cm_motion.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace constraint
{

  Cm_motion::Cm_motion(CAVIAR *fptr) : Constraint{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    velocity_steps = -1;
    angular_momentum_steps = -1;
    constraint_type = Constraint_t::Cm_motion;
  }

  Cm_motion::~Cm_motion() {}

  bool Cm_motion::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else if (string_cmp(t, "velocity_steps"))
      {
        GET_OR_CHOOSE_A_INT(velocity_steps, "", "")
        if (velocity_steps < 1)
          error->all(FC_FILE_LINE_FUNC_PARSE, "velocity_steps have to be non-negative.");
      }
      else if (string_cmp(t, "angular_momentum_steps"))
      {
        GET_OR_CHOOSE_A_INT(angular_momentum_steps, "", "")
        if (angular_momentum_steps < 1)
          error->all(FC_FILE_LINE_FUNC_PARSE, "angular_momentum_steps have to be non-negative.");
      }
      else
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
      }
    }
    return in_file;
  }

  void Cm_motion::fix_velocity(int64_t steps, bool &recalculate_temperature)
  { // step I
    // XXX there may be two cases. 1: all of particles, 2: a particle of a type

    FC_NULLPTR_CHECK(atom_data)

    if (steps % velocity_steps == 0)
    {
      recalculate_temperature = true;
      fix_linear_momentum();
    }

    if (steps % angular_momentum_steps == 0)
    {
      recalculate_temperature = true;
      fix_angular_momentum();
    }
  }

  void Cm_motion::fix_linear_momentum()
  {
    auto v_cm = atom_data->owned_velocity_cm();
    for (auto &&v : atom_data->atom_struct_owned.velocity)
    {
      v -= v_cm;
    }
  }

  void Cm_motion::fix_angular_momentum()
  {
    error->all(FC_FILE_LINE_FUNC, "not implemented.");
    // auto v_cm = atom_data->owned_angular_momentum_cm();
  }

} // constraint

CAVIAR_NAMESPACE_CLOSE
