
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

#include "caviar/objects/constraint/nose_hoover.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace constraint
{

  Nose_hoover::Nose_hoover(CAVIAR *fptr) : Constraint{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    dt = -1.0;
    settings_verified = false;
    kb = 1.0;
    mass = -1.0;
    tau = -1.0;
    type = 1;
    constraint_type = Constraint_t::Nose_hoover;
  }

  Nose_hoover::~Nose_hoover() {}

  bool Nose_hoover::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR
      if (string_cmp(t, "mass"))
      {
        GET_OR_CHOOSE_A_REAL(mass, "", "")
        if (mass <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "mass have to non-negative.");
      }
      else if (string_cmp(t, "dt"))
      {
        GET_OR_CHOOSE_A_REAL(dt, "", "")
        if (dt <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "dt have to non-negative.");
      }
      else if (string_cmp(t, "tau"))
      {
        GET_OR_CHOOSE_A_REAL(tau, "", "")
        if (tau <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "dt have to non-negative.");
      }
      else if (string_cmp(t, "type"))
      {
        GET_OR_CHOOSE_A_INT(type, "", "")
        if (type <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "type have to non-negative.");
      }
      else if (string_cmp(t, "kb"))
      {
        GET_OR_CHOOSE_A_REAL(kb, "", "")
        if (kb <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "kb have to non-negative.");
      }
      else if (string_cmp(t, "temperature"))
      {
        GET_OR_CHOOSE_A_REAL(temperature, "", "")
        if (temperature <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Temperature have to non-negative.");
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

  void Nose_hoover::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)

    if (dt < 0.0)
      error->all(FC_FILE_LINE_FUNC, "dt is not set.");

    switch (type)
    {
    case (1):
    {
      if (mass < 0.0)
        error->all(FC_FILE_LINE_FUNC, "mass is not set.");
    }
    break;

    case (2):
    {
      if (tau < 0.0)
        error->all(FC_FILE_LINE_FUNC, "tau is not set.");
    }
    break;

    default:
      error->all(FC_FILE_LINE_FUNC, "this Nose-Hoover type is not implemented.  Expected type 1 or 2");
    }
    settings_verified = true;
  }

  void Nose_hoover::apply_thermostat(int64_t timestep,  bool &recalculate_temperature)
  { // step I
    if (!settings_verified)
      verify_settings();
    recalculate_temperature = true;
    auto n_df = atom_data->degree_of_freedoms();

    // Nose formalism (real-time)
    auto g = n_df;

    // Nose-Hoover formalism (virtual-time)
    // auto g = n_df + 1;

    auto temp_inst = atom_data->temperature();

    if (temp_inst == 0)
    {
      output->warning("Temperature = 0. Nose_hoover thermostat step is ignored at "+ std::to_string (timestep));
    }

    switch (type)
    {
    case (1):
    {
      zeta_dot = (-kb * n_df * temp_inst / mass) * ((g * temperature / (n_df * temp_inst)) - 1.0);
    }
    break;

    case (2):
    {
      zeta_dot = (-1.0 / tau) * ((g * temperature / (n_df * temp_inst)) - 1.0);
    }
    break;

    default:
      error->all(FC_FILE_LINE_FUNC, "this Nose-Hoover type is not implemented. Expected type 1 or 2");
    }

    // a simple euler integration.
    zeta += dt * zeta_dot;

    // ------------------------------- // step II
    // We can put both of these function (step_I and step_II) into one if the Euler
    // integration is used.

    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &acc = atom_data->atom_struct_owned.acceleration;
    auto psize = acc.size();
    for (unsigned int i = 0; i < psize; ++i)
    {
      acc[i] += -zeta * vel[i];
    }
  }

} // constraint

CAVIAR_NAMESPACE_CLOSE
