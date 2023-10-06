
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

#include "caviar/objects/constraint/nve.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace constraint
{

  Nve::Nve(CAVIAR *fptr) : Constraint{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    energy_per_dof = -1.0;

    energy_tot = -1.0;

    temperature = -1.0;
    kb = -1.0;

    kbt = -1.0;
    constraint_type = Constraint_t::Nve;
  }

  Nve::~Nve() {}

  bool Nve::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "energy_per_dof"))
      {
        GET_OR_CHOOSE_A_REAL(energy_per_dof, "", "")
        if (energy_per_dof <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "energy_per_dof have to non-negative.");
      }
      else if (string_cmp(t, "energy_tot"))
      {
        GET_OR_CHOOSE_A_REAL(energy_tot, "", "")
        if (energy_tot <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "energy_tot have to non-negative.");
      }
      else if (string_cmp(t, "step"))
      {
        GET_OR_CHOOSE_A_INT(step, "", "")
        if (step <= 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "step have to non-negative.");
      }
      else if (string_cmp(t, "kb"))
      {
        GET_OR_CHOOSE_A_REAL(kb, "", "")
        if (kb <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "kb have to non-negative.");
      }
      else if (string_cmp(t, "kbt"))
      {
        GET_OR_CHOOSE_A_REAL(kbt, "", "")
        if (kbt <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "kbt have to non-negative.");
      }
      else if (string_cmp(t, "temperature"))
      {
        GET_OR_CHOOSE_A_REAL(temperature, "", "")
        if (temperature <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "temperature have to non-negative.");
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

  void Nve::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)

    bool kbt_set = (kbt > 0.0);

    bool temperature_set = (temperature > 0.0);
    bool kb_set = (kb > 0.0);
    bool kb_temperature_set = (temperature_set && kb_set);

    if ((temperature_set && kb_set) != (temperature_set || kb_set))
      error->all(FC_FILE_LINE_FUNC, "kb and temperature should be set together.");

    bool kbt_kb_temperature_set = (kb_temperature_set || kbt_set);

    bool energy_set = (energy_per_dof > 0.0 || energy_tot > 0.0);

    if (!(energy_set || kbt_kb_temperature_set))
      error->all(FC_FILE_LINE_FUNC, "energy or temperature is not set.");

    if (energy_set && kbt_kb_temperature_set)
      error->all(FC_FILE_LINE_FUNC, "Only one of 'energy' and 'temperature' can be set.");

    if (energy_set)
      if (energy_per_dof > 0.0 && energy_tot > 0.0)
        error->all(FC_FILE_LINE_FUNC, "Only one of 'energy_per_dof' and 'energy_tot' can be set.");

    if (kbt_kb_temperature_set)
    {
      if (kbt_set && kb_temperature_set)
        error->all(FC_FILE_LINE_FUNC, "Only one of 'kbt' and 'temperatue kb' can be set.");
    }
  }

  void Nve::apply_thermostat(int64_t timestep, bool &recalculate_temperature)
  {

    FC_OBJECT_VERIFY_SETTINGS

    if (timestep % step != 0) return;

    recalculate_temperature = true;

    auto energy_now = atom_data->kinetic_energy();

    double lambda = 1.0; // conversion coefficient

    if (energy_now == 0)
    {
      output->warning("energy = 0. NVE thermostat step is ignored at "+ std::to_string (timestep));
    }

    if (energy_per_dof > 0.0)
    {

      auto dof = atom_data->degree_of_freedoms();

      // energy_per_dof = lambda^2 * energy_now / dof
      lambda = std::sqrt(energy_per_dof * dof / energy_now);
    }
    else if (energy_tot > 0.0)
    {

      // energy_tot = lambda^2 * energy_now;
      lambda = std::sqrt(energy_tot / energy_now);
    }
    else if (kbt > 0.0)
    {

      auto dof = atom_data->degree_of_freedoms();

      // lambda^2 *energy_now = dof * kbt / 2

      lambda = std::sqrt(dof * kbt * 0.5 / energy_now);
    }
    else if (temperature > 0.0)
    {

      auto dof = atom_data->degree_of_freedoms();

      lambda = std::sqrt(dof * kb * temperature * 0.5 / energy_now);
    }

    auto &vel = atom_data->atom_struct_owned.velocity;

    for (auto &&v : vel)
      v *= lambda;
  }

} // constraint

CAVIAR_NAMESPACE_CLOSE
