
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Implementation of Shake algorithm is done by Ashkan Shahmoradi
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

#include "caviar/objects/constraint/shake.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"
#include "caviar/utility/interpreter_io_headers.h"

// #define ABS(x) (x<0 ? -x : x)

CAVIAR_NAMESPACE_OPEN

namespace constraint
{

  Shake::Shake(CAVIAR *fptr) : Constraint{fptr},
                               domain{nullptr},
                               dt{-1.0}, error_tolerance{1e-6},
                               initialized{false}
  {
    FC_OBJECT_INITIALIZE_INFO
    shake_type = 0;
    constraint_type = Constraint_t::Shake;
  }

  Shake::~Shake() {}

  bool Shake::read(caviar::interpreter::Parser *parser)
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
      else if (string_cmp(t, "shake_type"))
      {
        GET_OR_CHOOSE_A_INT(shake_type, "", "")
      }
      else if (string_cmp(t, "error_tolerance"))
      {
        GET_OR_CHOOSE_A_REAL(error_tolerance, "", "")
      }
      else if (string_cmp(t, "set_domain") || string_cmp(t, "domain"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }
    return in_file;
  }

  void Shake::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    atom_data->set_record_owned_position_old(true);

    FC_NULLPTR_CHECK(domain)

    if (dt <= 0.0)
      error->all(FC_FILE_LINE_FUNC, "dt have to be a positive number");
  }

  void Shake::apply_shake(int64_t)
  {

    FC_OBJECT_VERIFY_SETTINGS


    double virialConstraintLocal = 0;
    auto &pos = atom_data->atom_struct_owned.position;
    auto &atom_data_pos_old = atom_data->atom_struct_owned.position_old;

    double dt2 = dt * dt;

    //double test1 = 0;
    //double test2 = 0;
    double test3 = 0;


    for (unsigned int i = 0; i < atom_data->molecule_struct_owned.size(); i++)
    {

      if (atom_data->molecule_struct_owned[i].ghost)
        continue;

      auto &atomic_bond_vector = atom_data->molecule_struct_owned[i].atomic_bond_vector;

      auto Nc = atomic_bond_vector.size();
      if (Nc == 0)
        continue;

      std::vector<double> l(Nc, 0);

      double error_max = 1.0;
      int counter = -1;
      while (error_max > error_tolerance)
      {
        counter++;
        if (counter > iteration_max)
          error->all (FC_FILE_LINE_FUNC, "Iteration didn't converge in " + std::to_string(counter) + " iterations");

        for (unsigned int j = 0; j < atomic_bond_vector.size(); j++)
        {
          int id_1 = atomic_bond_vector[j].id_1, id_2 = atomic_bond_vector[j].id_2;
          int k1 = atom_data->atom_id_to_index[id_1], k2 = atom_data->atom_id_to_index[id_2];

          double d = atomic_bond_vector[j].length;

          double mass_inv_k1 = atom_data->atom_type_params.mass_inv[atom_data->atom_struct_owned.type[k1]];
          double mass_inv_k2 = atom_data->atom_type_params.mass_inv[atom_data->atom_struct_owned.type[k2]];

          auto dr = domain->fix_distance(pos[k1] - pos[k2]);

          auto r2 = dr * dr;

          auto dr_old = domain->fix_distance(atom_data_pos_old[k1] - atom_data_pos_old[k2]);

          auto dot = dr * dr_old;

          l[j] = (r2 - d * d) / (4 * dt2 * dot * (mass_inv_k1 + mass_inv_k2));

          auto f_coef = -2.0 * dt2 * l[j];

          auto fc = -f_coef * dr_old;

          auto force_val = -2.0 * l[j]* std::sqrt(dr_old*dr_old);

          //test1+= (force_val) * std::sqrt(r2);
          //test2+= (force_val) * std::sqrt(dot);
          test3+= (force_val) * std::sqrt(dr_old*dr_old);

          pos[k1] -= fc * mass_inv_k1;
          pos[k2] += fc * mass_inv_k2;
        }

        error_max = 0.0;
        for (unsigned int j = 0; j < atomic_bond_vector.size(); j++)
        {
          int k1 = atomic_bond_vector[j].id_1, k2 = atomic_bond_vector[j].id_2;

          auto d = atomic_bond_vector[j].length;

          auto d2 = d * d;

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
      // if (atom_data->get_pressure_process())
      // {
      //   for (unsigned int j = 0; j < atomic_bond_vector.size(); j++)
      //   {
      //     int id_1 = atomic_bond_vector[j].id_1, id_2 = atomic_bond_vector[j].id_2;
      //     int k1 = atom_data->atom_id_to_index[id_1], k2 = atom_data->atom_id_to_index[id_2];

      //     auto dr_old = domain->fix_distance(atom_data_pos_old[k1] - atom_data_pos_old[k2]);
      //     auto dr_old2 = dr_old * dr_old;
          
      //     virialConstraintLocal += -2.0 * l[j] * dr_old2; //
      //   }
      // }
    }
    //atom_data->virialConstraint += virialConstraintLocal;
    atom_data->virialConstraint += test3;






    // =====================
    //
    // =====================
    update_velocity_after_shake(-1);

  }

  void Shake::update_velocity_after_shake(int64_t)
  {
    // velocity_fix part
    // this fix has to be done only on the M-Shake molecules. If not, the normal
    // leap-frog step has to be enough.
    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &pos = atom_data->atom_struct_owned.position;
    auto &atom_data_pos_old = atom_data->atom_struct_owned.position_old;
    double dtinv = 1.0 / dt;

    for (unsigned int i = 0; i < atom_data->molecule_struct_owned.size(); i++)
    {

      if (atom_data->molecule_struct_owned[i].ghost)
        continue;

      auto &atomic_bond_vector = atom_data->molecule_struct_owned[i].atomic_bond_vector;

      for (unsigned int j = 0; j < atomic_bond_vector.size(); j++)
      { // XXX P.II
        auto id_1 = atomic_bond_vector[j].id_1;
        auto id_2 = atomic_bond_vector[j].id_2;
        int k1 = atom_data->atom_id_to_index[id_1], k2 = atom_data->atom_id_to_index[id_2];

        vel[k1] = domain->fix_distance(pos[k1] - atom_data_pos_old[k1]) * dtinv;
        vel[k2] = domain->fix_distance(pos[k2] - atom_data_pos_old[k2]) * dtinv;
      }
    }
  }

} // constraint

CAVIAR_NAMESPACE_CLOSE
