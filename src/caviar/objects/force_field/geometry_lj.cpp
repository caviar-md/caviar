
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

#include "caviar/objects/force_field/geometry_lj.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/shape.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/unique/time_function_3d.h"
#include <string>
#include <cmath>
#include <fstream>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Geometry_lj::Geometry_lj(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    wca = false;
    cutoff_list_activated = false;
    force_coef = 1.0;
  }

  Geometry_lj::~Geometry_lj()
  {
  }

  bool Geometry_lj::read(class caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "cutoff"))
      {
        GET_OR_CHOOSE_A_REAL(cutoff, "", "")
        if (cutoff < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Force field cutoff have to non-negative.");
      }
      else if (string_cmp(t, "cutoff_list"))
      {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT(cutoff_list)
        cutoff_list_activated = true;
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Epsilon have to be non-negative.");
      }
      else if (string_cmp(t, "wca"))
      {
        wca = true;
        cutoff_list_activated = true;
      }
      else if (string_cmp(t, "add_shape"))
      {
        FIND_OBJECT_BY_NAME(shape, it)
        shape.push_back(object_container->shape[it->second.index]);
        shape_type.push_back(0);
      }
      else if (string_cmp(t, "shape_type"))
      {
        unsigned int i = 0, t = 0;
        GET_OR_CHOOSE_A_INT(i, "", "")
        GET_OR_CHOOSE_A_INT(t, "", "")
        if (shape_type.size() < (unsigned int)i + 1)
          shape_type.resize(i + 1, 0);
        shape_type[i] = t;
      }
      else if (string_cmp(t, "force_coef"))
      {
        GET_OR_CHOOSE_A_REAL(force_coef, "", "")
      }
      else if (string_cmp(t, "epsilon_atom"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(epsilon_atom)
      }
      else if (string_cmp(t, "epsilon_wall"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(epsilon_wall)
      }
      else if (string_cmp(t, "sigma_atom"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(sigma_atom)
      }
      else if (string_cmp(t, "sigma_wall"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(sigma_wall)
      } /*else if (string_cmp(t,"epsilon")) {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT(epsilon)
        if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Epsilon have to be non-negative.");
      } else if (string_cmp(t,"sigma")) {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT(sigma)
        if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Sigma have to be non-negative.");
      } */
      else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else if (string_cmp(t, "set_position_offset"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        FC_CHECK_OBJECT_CLASS_NAME(unique, it, time_function_3d)
        unique::Time_function_3d *a = dynamic_cast<unique::Time_function_3d *>(object_container->unique[it->second.index]);
        position_offset = a;
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Geometry_lj::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)

    auto shape_size = shape.size();
    if (shape_size == 0)
    {
      output->warning("Geometry_lj:: shape.size()==0");
    }

    auto epsilon_wall_size = epsilon_wall.size();
    auto sigma_wall_size = sigma_wall.size();
    auto wall_max_size = (epsilon_wall_size > sigma_wall_size ? epsilon_wall_size : sigma_wall_size);
    if (epsilon_wall_size != sigma_wall_size)
    {
      output->warning("Geometry_lj:: (epsilon_wall_size != sigma_wall_size)");
      if (epsilon_wall_size != wall_max_size)
      {
        epsilon_wall.resize(wall_max_size, 0);
      }
      else
      {
        sigma_wall.resize(wall_max_size, 0);
      }
    }

    unsigned shape_type_max = 0;
    for (auto &&t : shape_type)
    {
      if ((int)shape_type_max < t)
        shape_type_max = t;
    }

    if (wall_max_size < shape_type_max)
    {
      output->warning("Geometry_lj:: (wall_max_size < shape_type_max)");
      epsilon_wall.resize(shape_type_max, 0);
      sigma_wall.resize(shape_type_max, 0);
      wall_max_size = shape_type_max;
    }

    auto epsilon_atom_size = epsilon_atom.size();
    auto sigma_atom_size = sigma_atom.size();
    auto atom_max_size = (epsilon_atom_size > sigma_atom_size ? epsilon_atom_size : sigma_atom_size);
    if (epsilon_atom_size != sigma_atom_size)
    {
      output->warning("Geometry_lj:: (epsilon_atom_size != sigma_atom_size)");
      if (epsilon_atom_size != atom_max_size)
      {
        epsilon_atom.resize(atom_max_size, 0);
      }
      else
      {
        sigma_atom.resize(atom_max_size, 0);
      }
    }

    unsigned atom_data_type_max = 0;
    for (auto &&t : atom_data->atom_struct_owned.type)
    {
      if (atom_data_type_max < t)
        atom_data_type_max = t;
    }

    if (atom_max_size < atom_data_type_max + 1)
    {
      output->warning("Geometry_lj:: (atom_max_size < atom_data_type_max + 1)");
      epsilon_atom.resize(atom_data_type_max + 1, 0);
      sigma_atom.resize(atom_data_type_max + 1, 0);
      atom_max_size = atom_data_type_max + 1;
    }

    epsilon.resize(wall_max_size);
    sigma.resize(wall_max_size);

    for (auto &&e : epsilon)
      e.resize(atom_max_size, 0);

    for (auto &&s : sigma)
      s.resize(atom_max_size, 0);

    if (true)
    {
      for (unsigned int i = 0; i < wall_max_size; ++i)
      {
        for (unsigned int j = 0; j < atom_max_size; ++j)
        {
          auto sigma_ij = 0.5 * (sigma_wall[i] + sigma_atom[j]);
          auto epsilon_ij = std::sqrt(epsilon_wall[i] * epsilon_atom[j]);
          sigma[i][j] = sigma_ij;
          epsilon[i][j] = epsilon_ij;
        }
      }
    }

    // Week-Chandler-Anderson (WCA) potential activated.
    if (wca)
    {
      auto cut_coef = std::pow(2.0, 1.0 / 6.0); // ONLY for LJ 6-12

      cutoff_list.resize(wall_max_size);
      for (auto &&c : cutoff_list)
        c.resize(atom_max_size, 0);

      for (unsigned int i = 0; i < wall_max_size; ++i)
      {
        auto type_s = shape_type[i];
        for (unsigned int j = 0; j < atom_max_size; ++j)
        {
          auto cut = cut_coef * sigma[type_s][j];
          std::cout << cut << "\n";
          cutoff_list[type_s][j] = cut;
        }
      }
    }
    else if (cutoff_list_activated)
    {
      if (cutoff_list.size() < shape_type_max)
      {
        cutoff_list.resize(shape_type_max);
        output->warning("Geometry_lj:: (cutoff_list.size() < shape_type_max)");
      }
      for (auto &&c : cutoff_list)
      {
        if (c.size() < atom_max_size)
        {
          c.resize(atom_max_size, 0);
          output->warning("Geometry_lj:: (cc.size() < atom_max_size)");
        }
      }
    }
  }

  void Geometry_lj::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    Vector<double> p_o{0, 0, 0};
    if (position_offset != nullptr)
      p_o = position_offset->current_value;

    const auto &pos = atom_data->atom_struct_owned.position;
    auto &acc = atom_data->atom_struct_owned.acceleration;

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < pos.size(); ++i)
    {
      const auto type_i = atom_data->atom_struct_owned.type[i];
      const auto mass_inv_i = atom_data->atom_type_params.mass_inv[type_i];

      for (unsigned int j = 0; j < shape.size(); ++j)
      {
        auto type_s = shape_type[j];

        Vector<Real_t> contact_vector{0, 0, 0};

        double c = cutoff;

        if (cutoff_list_activated)
          c = cutoff_list[type_s][type_i];

        bool is_in_contact = shape[j]->in_contact(pos[i] - p_o, c, contact_vector);

        // if distance to the wall is less than cutoff.
        if (is_in_contact)
        {

          auto contact_vector_sq = contact_vector * contact_vector;

          // TODO not efficient. Fix it. maybe import it from in_contact()
          auto contact_vector_sqrt = std::sqrt(contact_vector_sq);

          auto contact_vector_unit = contact_vector / contact_vector_sqrt;
          auto d = c - contact_vector_sqrt; // distance between wall's atom and atom.

          auto dr = contact_vector_unit * d;
          auto dr_sq = dr * dr;

          auto eps_ij = epsilon[type_s][type_i];
          auto sigma_ij = sigma[type_s][type_i];

          /*  //  XXX lj 10-4

                  auto eps_ij = epsilon [0][type_i];
                  auto sig_ij = sigma   [0][type_i];
                  auto sig_ij_inv = 1.0 / sig_ij;
                  auto force_coef = 2.0 * 3.141592653589 * eps_ij * sig_ij * sig_ij * force_coef;

                  auto dr_sq     = contact_vector * contact_vector;
          //      auto dr_sq_inv = 1.0 / dr_sq;
                  auto r_norm   = std::sqrt (dr_sq);
                  auto r_inv    = 1.0 / r_norm;

                  if (dr_sq > c*c) continue;

                  auto rho      = sig_ij * r_inv;
                  auto rho_sq   = rho * rho;
                  auto rho_5    = rho_sq * rho_sq * rho;

                  auto r_cut        = c;
          //      auto r_cut_sq     = r_cut * r_cut;
          //      auto r_cut_sq_inv = 1.0 / r_cut_sq;
                  auto r_cut_norm   = r_cut;
                  auto r_cut_inv    = 1.0 / r_cut_norm;

                  auto rho_cut      = sig_ij * r_cut_inv;
                  auto rho_cut_inv  = 1.0 / rho_cut;
                  auto rho_cut_sq   = rho_cut * rho_cut;
                  auto rho_cut_5    = rho_cut_sq * rho_cut_sq * rho_cut;

                  auto force_norm     = force_coef * 4.0 *  ( - (rho_5*rho_5*r_inv)               + (rho_5*sig_ij_inv));
                  auto force_norm_cut = force_coef * 4.0 *  ( - (rho_cut_5*rho_cut_5*rho_cut_inv) + (rho_cut_5*sig_ij_inv));

                  auto force = (force_norm - force_norm_cut) * r_inv * contact_vector;

                  acc[i] += force*mass_inv_i;
          */

          // XXX lj 6-12

          auto r_c_sq_inv = 1 / (c * c);
          auto rho_c_sq_inv = sigma_ij * sigma_ij * r_c_sq_inv;
          auto rho_c_6_inv = rho_c_sq_inv * rho_c_sq_inv * rho_c_sq_inv;
          auto rho_c_12_inv = rho_c_6_inv * rho_c_6_inv;

          auto dr_sq_inv = 1 / dr_sq;
          auto rho_sq_inv = sigma_ij * sigma_ij * dr_sq_inv;
          auto rho_6_inv = rho_sq_inv * rho_sq_inv * rho_sq_inv;
          auto rho_12_inv = rho_6_inv * rho_6_inv;

          auto force = force_coef * 4 * eps_ij * (-12 * rho_12_inv * dr_sq_inv + 6 * rho_6_inv * dr_sq_inv + +12 * rho_c_12_inv * r_c_sq_inv - 6 * rho_c_6_inv * r_c_sq_inv) * dr;
          // std::cout << "f : " << force << "\n";

          acc[i] += force * mass_inv_i;
        }
      }
    }
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
