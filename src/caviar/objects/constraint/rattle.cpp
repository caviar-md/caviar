
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Implementation of Rattle algorithm is done by Ashkan Shahmoradi
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

#include "caviar/objects/constraint/rattle.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace constraint
{

  Rattle::Rattle(CAVIAR *fptr) : Constraint{fptr},
                                 domain{nullptr},
                                 dt{-1.0}, error_tolerance{1e-6},
                                 domain_dh{Vector<double>{0, 0, 0}},
                                 domain_bc{Vector<int>{0, 0, 0}},
                                 initialized{false}
  {
    FC_OBJECT_INITIALIZE_INFO
    constraint_type = Constraint_t::Rattle;
  }

  Rattle::~Rattle() {}

  bool Rattle::read(caviar::interpreter::Parser *parser)
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
      else if (string_cmp(t, "dt"))
      {
        GET_OR_CHOOSE_A_REAL(dt, "", "")
      }
      else if (string_cmp(t, "iteration_max"))
      {
        GET_OR_CHOOSE_A_INT(iteration_max, "", "")
      }
      else if (string_cmp(t, "error_tolerance"))
      {
        GET_OR_CHOOSE_A_REAL(error_tolerance, "", "")
      }
      else if (string_cmp(t, "set_domain") || string_cmp(t, "domain"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
        domain_dh = 0.5 * (domain->upper_global - domain->lower_global);
        domain_bc = domain->boundary_condition;
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }
    return in_file;
  }

  void Rattle::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    atom_data->set_record_owned_position_old(true);
    FC_NULLPTR_CHECK(domain)
    if (dt <= 0.0)
      error->all(FC_FILE_LINE_FUNC, "dt have to be a positive number");
  }

  void Rattle::apply_shake(int64_t)
  {
    FC_OBJECT_VERIFY_SETTINGS

    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &pos = atom_data->atom_struct_owned.position;
    auto &atom_data_pos_old = atom_data->atom_struct_owned.position_old;

    for (unsigned int i = 0; i < atom_data->molecule_struct_owned.size(); i++)
    {

      if (atom_data->molecule_struct_owned[i].ghost)
        continue;

      auto &atomic_bond_vector = atom_data->molecule_struct_owned[i].atomic_bond_vector;

      auto Nc = atomic_bond_vector.size();
      if (Nc == 0)
        continue;
      //std::vector<double> C(Nc, 0);

      double error_max = 1.0;
      int counter = -1;

      while (error_max > error_tolerance)
      {

        counter++;
        if (counter > iteration_max)
          error->all (FC_FILE_LINE_FUNC, "Iteration didn't converge in " + std::to_string(counter) + " iterations");

        for (unsigned int j = 0; j < atomic_bond_vector.size(); j++)
        {


          auto id_1 = atomic_bond_vector[j].id_1;
          auto id_2 = atomic_bond_vector[j].id_2;
          int k1 = atom_data->atom_id_to_index[id_1];
          int k2 = atom_data->atom_id_to_index[id_2];

          auto d = atomic_bond_vector[j].length;

          auto mass_inv_k1 = atom_data->atom_type_params.mass_inv[atom_data->atom_struct_owned.type[k1]];
          auto mass_inv_k2 = atom_data->atom_type_params.mass_inv[atom_data->atom_struct_owned.type[k2]];

          auto dr = domain->fix_distance(pos[k1] - pos[k2]);

          auto dr_old = domain->fix_distance(atom_data_pos_old[k1] - atom_data_pos_old[k2]);

          auto lambda = (dr * dr - (d * d)) / (2.0 * (mass_inv_k1 + mass_inv_k2) * (dr * dr_old));

          pos[k1] -= mass_inv_k1 * lambda * dr_old;

          //dr_old = domain->fix_distance(atom_data_pos_old[k2] - atom_data_pos_old[k1]); // Note (k2 - k1)
          dr_old = - dr_old; // ??????????? Check the paper

          pos[k2] -= mass_inv_k2 * lambda * dr_old;

          dr = domain->fix_distance(pos[k1] - pos[k2]);

          auto dv = vel[k1] - vel[k2];

          auto etha = (dr * dv) / ((mass_inv_k1 + mass_inv_k2) * d * d);

          vel[k1] -= mass_inv_k1 * etha * dr;

          vel[k2] += mass_inv_k2 * etha * dr; // Note the Plus sign
        }

        error_max = 0.0;
        for (unsigned int j = 0; j < atomic_bond_vector.size(); j++)
        {
          int k1 = atomic_bond_vector[j].id_1, k2 = atomic_bond_vector[j].id_2;

          auto d = atomic_bond_vector[j].length;

          auto d2 = d*d;

          auto dr = domain->fix_distance(pos[k1] - pos[k2]);

          auto r2 = dr * dr;

          double err = (r2 - d2) / (2 * d2);

          double abs_err = (err < 0) ? -err : err;

          if (error_max < abs_err)
            error_max = abs_err;
        }

      }
      // =====================
      // Virial Calculations
      // =====================
      if (atom_data->get_pressure_process())
      {
        error->all(FC_FILE_LINE_FUNC,"Pressure (Virial) calculation is Not implemented for this constraint");
      }
    }
  }

} // constraint

CAVIAR_NAMESPACE_CLOSE
