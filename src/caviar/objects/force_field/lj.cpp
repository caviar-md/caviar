
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

#include "caviar/objects/force_field/lj.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"
#include <cmath>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Lj::Lj(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    jump_fix = false;
    monitor_jump = false;
    jump_tol = 1e-6;
    wca = false;
    make_off_diagonal_vectors = false;
    cutoff_list_activated = false;
    input_by_array = false;
    input_by_atom = false;
  }

  bool Lj::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR
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
      else if (string_cmp(t, "epsilon"))
      {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT(epsilon)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Epsilon have to be non-negative.");
        input_by_array = true;
      }
      else if (string_cmp(t, "sigma"))
      {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT(sigma)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Sigma have to be non-negative.");
        input_by_array = true;
      }
      else if (string_cmp(t, "lambda_s"))
      {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT_R(lambda_s, 1.0)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Sigma have to be non-negative.");
        lambda_s_is_set = true;
      }
      else if (string_cmp(t, "lambda_e"))
      {
        GET_A_STDVECTOR_STDVECTOR_REAL_ELEMENT_R(lambda_e, 1.0)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Sigma have to be non-negative.");
        lambda_e_is_set = true;
      }
      else if (string_cmp(t, "epsilon_atom"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(epsilon_atom)
        input_by_atom = true;
      }
      else if (string_cmp(t, "sigma_atom"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(sigma_atom)
        input_by_atom = true;
      }
      else if (string_cmp(t, "set_neighborlist") || string_cmp(t, "neighborlist"))
      {
        FIND_OBJECT_BY_NAME(neighborlist, it)
        neighborlist = object_container->neighborlist[it->second.index];
      }
      else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else if (string_cmp(t, "set_domain") || string_cmp(t, "domain"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
      }
      else if (string_cmp(t, "jump_tol"))
      {
        // can be used only if there's a random force. It should happen  rarely.
        // it fixes the enormous amount of LJ to a lesser value.
        jump_fix = true;
        GET_OR_CHOOSE_A_REAL(jump_tol, "", "")
      }
      else if (string_cmp(t, "monitor_jump"))
      {
        monitor_jump = true;
      }
      else if (string_cmp(t, "ignore_intra_molecule"))
      {
        ignore_intra_molecule = true;
      }
      else if (string_cmp(t, "wca"))
      {
        wca = true;
        cutoff_list_activated = true;
      }
      else if (string_cmp(t, "make_off_diagonal_vectors"))
      {
        make_off_diagonal_vectors = true;
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Lj::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    // FC_NULLPTR_CHECK(domain)
    FC_NULLPTR_CHECK(neighborlist)
    my_mpi_rank = atom_data->get_mpi_rank();
    if (input_by_atom && input_by_array)
    {
      error->all(FC_FILE_LINE_FUNC, "(input_by_atom && input_by_array)");
    }

    unsigned max_size = 0;

    unsigned atom_data_type_max = 0;
    for (auto &&t : atom_data->atom_struct_owned.type)
    {
      // std::cout << "t: " << t << std::endl;
      if (atom_data_type_max < t)
        atom_data_type_max = t;
    }

    if (input_by_atom)
    {

      auto epsilon_atom_size = epsilon_atom.size();
      auto sigma_atom_size = sigma_atom.size();
      auto atom_max_size = (epsilon_atom_size > sigma_atom_size ? epsilon_atom_size : sigma_atom_size);
      if (epsilon_atom_size != sigma_atom_size)
      {
        output->warning("lj:: (epsilon_atom_size != sigma_atom_size)");
        if (epsilon_atom_size != atom_max_size)
        {
          epsilon_atom.resize(atom_max_size, 0);
        }
        else
        {
          sigma_atom.resize(atom_max_size, 0);
        }
      }

      if (atom_max_size < atom_data_type_max + 1)
      {
        output->warning("lj:: (atom_max_size < atom_data_type_max + 1)"); // if the input is done by epsilon rather than epsilon_atom,
                                                                          // we will have this warning.
        epsilon_atom.resize(atom_data_type_max + 1, 0);
        sigma_atom.resize(atom_data_type_max + 1, 0);
        atom_max_size = atom_data_type_max + 1;
      }

      epsilon.resize(atom_max_size);
      sigma.resize(atom_max_size);

      for (auto &&e : epsilon)
        e.resize(atom_max_size, 0);

      for (auto &&s : sigma)
        s.resize(atom_max_size, 0);

      for (unsigned int i = 0; i < atom_max_size; ++i)
      {
        for (unsigned int j = 0; j < atom_max_size; ++j)
        {
          if (i == j)
          {
            sigma[i][j] = sigma_atom[i];
            epsilon[i][j] = epsilon_atom[i];
            continue;
          }
          auto sigma_ij = 0.5 * (sigma_atom[i] + sigma_atom[j]);
          auto epsilon_ij = std::sqrt(epsilon_atom[i] * epsilon_atom[j]);
          sigma[i][j] = sigma_ij;
          epsilon[i][j] = epsilon_ij;
        }
      }

      max_size = atom_max_size;
    }
    else if (input_by_array)
    {

      auto epsilon_size = epsilon.size();
      auto sigma_size = sigma.size();
      auto es_max_size = (epsilon_size > sigma_size ? epsilon_size : sigma_size);

      if (epsilon_size != sigma_size)
      {
        output->warning("lj:: (epsilon_size != sigma_size)");
        if (epsilon_size != es_max_size)
        {
          epsilon.resize(es_max_size);
        }
        else
        {
          sigma.resize(es_max_size);
        }
      }

      if (es_max_size < atom_data_type_max + 1)
      {
        output->warning("lj:: (es_max_size < atom_data_type_max + 1)");

        es_max_size = atom_data_type_max + 1;

        epsilon.resize(es_max_size);
        sigma.resize(es_max_size);
      }

      for (auto &&e : epsilon)
      {
        if (e.size() != es_max_size)
        {
          output->warning("lj:: (e.size != es_max_size)");
          e.resize(es_max_size, 0);
        }
      }

      for (auto &&s : sigma)
      {
        if (s.size() != es_max_size)
        {
          output->warning("s.size != es_max_size");
          s.resize(es_max_size, 0);
        }
      }

      max_size = es_max_size;
    }

    // Week-Chandler-Anderson (WCA) potential activated.
    if (wca)
    {
      auto cut_coef = std::pow(2.0, 1.0 / 6.0); // ONLY for LJ 6-12

      cutoff_list.resize(max_size);
      for (auto &&c : cutoff_list)
        c.resize(max_size, 0);

      for (unsigned int i = 0; i < max_size; ++i)
      {
        for (unsigned int j = 0; j < max_size; ++j)
        {
          auto cut = cut_coef * sigma[i][j];
          std::cout << "[INF] lj wca cutoff " << i << "," << j << " : " << cut << "\n";
          cutoff_list[i][j] = cut;
        }
      }
    }
    else if (cutoff_list_activated)
    {
      if (cutoff_list.size() < max_size)
      {
        cutoff_list.resize(max_size);
        output->warning("lj:: (cutoff_list.size() < atom_max_size)");
      }
      for (auto &&c : cutoff_list)
      {
        if (c.size() < max_size)
        {
          c.resize(max_size, 0);
          output->warning("lj:: (cc.size() < atom_max_size)");
        }
      }
    }

    if (lambda_s_is_set)
    {
      if (lambda_s.size() <= atom_data_type_max) lambda_s.resize(atom_data_type_max+1);
      for (unsigned int i = 0; i < lambda_s.size(); ++i)
      {
        if (lambda_s[i].size() <= atom_data_type_max) lambda_s[i].resize(atom_data_type_max+1, 1.0);
      }

      for (unsigned int i = 0; i < lambda_s.size(); ++i)
      {
        for (unsigned int j = 0; j < lambda_s[i].size(); ++j)
        {
          sigma[i][j] *= lambda_s[i][j];
        }
      }
    }
    if (lambda_e_is_set)
    {
      if (lambda_e.size() <= atom_data_type_max) lambda_e.resize(atom_data_type_max+1);
      for (unsigned int i = 0; i < lambda_e.size(); ++i)
      {
        if (lambda_e[i].size() <= atom_data_type_max) lambda_e[i].resize(atom_data_type_max+1, 1.0);
      }

      for (unsigned int i = 0; i < lambda_e.size(); ++i)
      {
        for (unsigned int j = 0; j < lambda_e[i].size(); ++j)
        {
          epsilon[i][j] *= lambda_e[i][j];
        }
      }
    }
    /*
    for (unsigned int i = 0; i < max_size; ++i) {
    for (unsigned int j = 0; j < max_size; ++j) {
      std::cout << "e " << i << " , " << j << " : " << epsilon[i][j] << std::endl;
      std::cout << "s " << i << " , " << j << " : " << sigma[i][j] << std::endl;
    }
    }
    */
  }

  void Lj::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS


    auto cutoff_sq = cutoff * cutoff;
    auto c = cutoff_sq;
    const auto &nlist = neighborlist->neighlist;

    double virialLocal = 0;

    // auto &mol_index = atom_data->atom_struct_owned.molecule_index;

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < nlist.size(); ++i)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif

      const auto &pos_i = atom_data->atom_struct_owned.position[i];
      const auto type_i = atom_data->atom_struct_owned.type[i];
      const auto mass_inv_i = atom_data->atom_type_params.mass_inv[type_i];
      //const auto mol_index_i = mol_index[i];
      int id_i = atom_data->atom_struct_owned.id[i];
      const auto mol_index_i = atom_data->atom_struct_owned.molecule_index[i];

      for (auto j : nlist[i])
      {
        bool is_ghost = j >= nlist.size();
        Vector<Real_t> pos_j;
        Real_t type_j, mass_inv_j;
        int id_j;
        int mol_index_j;

        if (is_ghost)
        {
          j -= nlist.size();
          pos_j = atom_data->atom_struct_ghost.position[j];
          type_j = atom_data->atom_struct_ghost.type[j];
          id_j = atom_data->atom_struct_ghost.id[j];
          mol_index_j = atom_data->atom_struct_ghost.molecule_index[j];
        }
        else
        {
          pos_j = atom_data->atom_struct_owned.position[j];
          type_j = atom_data->atom_struct_owned.type[j];
          id_j = atom_data->atom_struct_owned.id[j];
          mol_index_j = atom_data->atom_struct_owned.molecule_index[j];
        }


        if (ignore_intra_molecule && (mol_index_i == mol_index_j))        
          continue;        

        const auto dr = pos_j - pos_i;

        auto dr_sq = dr * dr;

        if (cutoff_list_activated)
          c = cutoff_list[type_i][type_j] * cutoff_list[type_i][type_j];

        if (dr_sq > c)
          continue;

        if (jump_fix)
        {

          if (monitor_jump)
            if (dr_sq < jump_tol)
              std::cout << "\ndr_sq: " << dr_sq << " ";

          if (dr_sq < jump_tol)
            dr_sq = jump_tol;
        }

        const auto eps_ij = epsilon[type_i][type_j];
        if (eps_ij == 0)
          continue; // used in water simulations. Hydrogen in water models usiually has zero value.
        const auto sigma_ij = sigma[type_i][type_j];

        auto r_c_sq_inv = 1 / c;
        auto rho_c_sq_inv = sigma_ij * sigma_ij * r_c_sq_inv;
        auto rho_c_6_inv = rho_c_sq_inv * rho_c_sq_inv * rho_c_sq_inv;
        auto rho_c_12_inv = rho_c_6_inv * rho_c_6_inv;
        auto dr_sq_inv = 1 / dr_sq;
        auto rho_sq_inv = sigma_ij * sigma_ij * dr_sq_inv;
        auto rho_6_inv = rho_sq_inv * rho_sq_inv * rho_sq_inv;
        auto rho_12_inv = rho_6_inv * rho_6_inv;

        auto forceCoef = 4 * eps_ij * (-12 * rho_12_inv * dr_sq_inv + 6 * rho_6_inv * dr_sq_inv + 12 * rho_c_12_inv * r_c_sq_inv - 6 * rho_c_6_inv * r_c_sq_inv);
        //if (i < j)
        //if ((mol_index_i == -1) || (mol_index_i != mol_index_j))
        if (id_i < id_j)
        {
          virialLocal +=  -forceCoef * dr_sq;
        }
        auto force = forceCoef * dr;


        atom_data->atom_struct_owned.acceleration[i] += force * mass_inv_i;
        
        // if (atom_data->pressure_process)
        //   atom_data->add_to_pressure(force*dr);

        if (!is_ghost)
        {
          mass_inv_j = atom_data->atom_type_params.mass_inv[type_j];

#ifdef CAVIAR_WITH_OPENMP
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[j].x -= force.x * mass_inv_j;
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[j].y -= force.y * mass_inv_j;
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[j].z -= force.z * mass_inv_j;
#else
          atom_data->atom_struct_owned.acceleration[j] -= force * mass_inv_j;
#endif
        }
      }
    }
    //std::cout << "lj virialLocal " << virialLocal << std::endl;
    atom_data->virialForce += virialLocal;

  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
