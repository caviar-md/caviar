
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

#include "caviar/objects/constraint/berendsen.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace constraint
{

  Berendsen::Berendsen(CAVIAR *fptr) : Constraint{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    coupling = -1.0;
    dt = -1.0;
    constraint_type = Constraint_t::Berendsen;
  }

  Berendsen::~Berendsen() {}

  bool Berendsen::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "temperature"))
      {
        GET_OR_CHOOSE_A_REAL(temperature, "", "")
        if (temperature <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Temperature have to non-negative.");
      }
      else if (string_cmp(t, "coupling"))
      {
        GET_OR_CHOOSE_A_REAL(coupling, "", "")
        if (coupling <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "coupling have to non-negative.");
      }
      else if (string_cmp(t, "dt"))
      {
        GET_OR_CHOOSE_A_REAL(dt, "", "")
        if (dt <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "dt have to non-negative.");
      }
      else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
      }
    }
    return in_file;
  }

  void Berendsen::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)

    if (dt < 0.0)
      error->all(FC_FILE_LINE_FUNC, "dt is not set.");
    if (coupling < 0.0)
      error->all(FC_FILE_LINE_FUNC, "coupling is not set.");
  }

  void Berendsen::apply_on_velocity(int64_t)
  { // step I

    FC_OBJECT_VERIFY_SETTINGS

    auto t = atom_data->temperature();

    auto &vel = atom_data->atom_struct_owned.velocity;

    double lambda = std::sqrt(1.0 + (dt / coupling) * ((temperature / t) - 1.0));

    for (auto &&v : vel)
      v *= lambda;
  }

} // constraint

CAVIAR_NAMESPACE_CLOSE
