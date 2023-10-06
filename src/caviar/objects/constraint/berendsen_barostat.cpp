
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

#include "caviar/objects/constraint/berendsen_barostat.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/force_field.h"
#include "caviar/objects/domain.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace constraint
{

  Berendsen_barostat::Berendsen_barostat(CAVIAR *fptr) : Constraint{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    tp = -1.0;
    dt = -1.0;
    constraint_type = Constraint_t::Berendsen;
    domain = nullptr;
  }

  Berendsen_barostat::~Berendsen_barostat() {}

  bool Berendsen_barostat::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "tp"))
      {
        GET_OR_CHOOSE_A_REAL(tp, "", "")
        if (tp <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "tp have to non-negative.");
      }
      else if (string_cmp(t, "kappa"))
      {
        GET_OR_CHOOSE_A_REAL(kappa, "", "")
        if (kappa <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "kappa have to non-negative.");
      }
      else if (string_cmp(t, "pressure"))
      {
        GET_OR_CHOOSE_A_REAL(pressure, "", "")
        if (pressure <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "pressure have to non-negative.");
      }
      else if (string_cmp(t, "dt"))
      {
        GET_OR_CHOOSE_A_REAL(dt, "", "")
        if (dt <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "dt have to non-negative.");
      }
      else if (string_cmp(t, "step"))
      {
        GET_OR_CHOOSE_A_INT(step, "", "")
        if (step <= 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "step have to non-negative.");
      }
      else if (string_cmp(t, "xi_max"))
      {
        GET_OR_CHOOSE_A_REAL(xi_max, "", "")
        if (xi_max <= 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "xi_max have to non-negative.");
      }
      else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else if (string_cmp(t, "add_force_field") || string_cmp(t, "force_field"))
      {
        FIND_OBJECT_BY_NAME(force_field, it)
        force_field.push_back(object_container->force_field[it->second.index]);
      }
      else if (string_cmp(t, "set_domain") || string_cmp(t, "domain"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
        scale_axis = domain->boundary_condition;
      }
        else if (string_cmp(t, "scale_axis") || string_cmp(t, "bc"))
      {
        GET_OR_CHOOSE_A_INT_3D_VECTOR(scale_axis, "", "")
      }
      else
      {
        error->all(FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
      }
    }
    return in_file;
  }

  void Berendsen_barostat::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(domain)

    if (dt < 0.0)
      error->all(FC_FILE_LINE_FUNC, "dt is not set.");
    if (tp <= 0.0)
      error->all(FC_FILE_LINE_FUNC, "tp is not set.");
    if (kappa < 0.0)
      error->all(FC_FILE_LINE_FUNC, "kappa is not set.");
    if (pressure < 0.0)
      error->all(FC_FILE_LINE_FUNC, "pressure is not set.");
  }

  void Berendsen_barostat::apply_barostat(int64_t timestep, bool &fix_position_needed)
  {     
    if (timestep % step != 0) return;
    FC_OBJECT_VERIFY_SETTINGS

    auto p = atom_data->pressure();
    //std::cout << "p: " << p << std::endl;

    double xi = std::pow( 1.0 - kappa * (dt / tp) * (pressure - p), 0.333333333333333 );
    
    //std::cout << "xi_calc: " << xi << std::endl;

    double xi_low = 1 - xi_max;
    double xi_high = 1 + xi_max;

    // to avoid crashes due to large or small scalings
    if (xi < xi_low) xi = xi_low ;
    if (xi > xi_high) xi = xi_high;

    //std::cout << "xi_fix: " << xi << std::endl;

    domain-> scale_position(xi, scale_axis);

    atom_data->scale_position(xi, scale_axis);

    for (unsigned int i = 0; i < force_field.size(); ++i)
      force_field[i]->scale_position(xi, scale_axis);
    
    fix_position_needed = true;

  }

} // constraint

CAVIAR_NAMESPACE_CLOSE
