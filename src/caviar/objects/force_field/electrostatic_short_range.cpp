
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

#include "caviar/objects/force_field/electrostatic_short_range.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Electrostatic_short_range::Electrostatic_short_range(CAVIAR *fptr) : Force_field{fptr},
                                                                       k_electrostatic{1.0}, external_field{Vector<double>{0, 0, 0}},
                                                                       beta{1.0}, initialized{false}
  {
    FC_OBJECT_INITIALIZE_INFO
    cutoff = -1.0;
  }

  bool Electrostatic_short_range::read(caviar::interpreter::Parser *parser)
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
      else if (string_cmp(t, "external_field"))
      {
        GET_OR_CHOOSE_A_REAL_3D_VECTOR(external_field, "", "");
        // if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Epsilon have to be non-negative.");
      }
      else if (string_cmp(t, "k_electrostatic"))
      {
        GET_OR_CHOOSE_A_REAL(k_electrostatic, "", "")
        if (k_electrostatic < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "k_electrostatic has to be non-negative.");
      }
      else if (string_cmp(t, "beta"))
      {
        GET_OR_CHOOSE_A_REAL(beta, "", "")
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
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Electrostatic_short_range::initialize()
  {
    initialized = true;
    C = 1.0 / (beta * std::pow(cutoff, beta + 1.0));
    D = -(beta + 1.0) / (beta * cutoff);
    cutoff_sq = cutoff * cutoff;
  }

  void Electrostatic_short_range::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    // FC_NULLPTR_CHECK(domain)
    FC_NULLPTR_CHECK(neighborlist)
    if (cutoff < 0.0)
      error->all(FC_FILE_LINE_FUNC, "Force field cutoff have to non-negative.");
    my_mpi_rank = atom_data->get_mpi_rank();
  }

  void Electrostatic_short_range::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS
    if (!initialized)
      initialize();

    // auto p18 = atom_data->atom_struct_owned.position[18];
    // auto p19 = atom_data->atom_struct_owned.position[19];
    // auto p20 = atom_data->atom_struct_owned.position[20];
        // std::cout << "==================== " << std::endl;

    // std::cout << "p18:" << p18.x << " " << p18.y << " " << p18.z << std::endl;
    // std::cout << "p19:" << p19.x << " " << p19.y << " " << p19.z << std::endl;
    // std::cout << "p20:" << p20.x << " " << p20.y << " " << p20.z << std::endl;
        // std::cout <<"--------------------- " << std::endl;

    for (int i = 0; i < atom_data->atom_struct_ghost.position.size(); ++i)
    {
      if (atom_data->atom_struct_ghost.id[i] == 18)
      {
        auto g18 = atom_data->atom_struct_ghost.position[i];
        // std::cout << "g18:" << g18.x << " " << g18.y << " " << g18.z << std::endl;
      }
      if (atom_data->atom_struct_ghost.id[i] == 19)
      {
        auto g19 = atom_data->atom_struct_ghost.position[i];
        // std::cout << "g19:" << g19.x << " " << g19.y << " " << g19.z << std::endl;
      }
      if (atom_data->atom_struct_ghost.id[i] == 20)
      {
        auto g20 = atom_data->atom_struct_ghost.position[i];
        // std::cout << "g20:" << g20.x << " " << g20.y << " " << g20.z << std::endl;
      }
    }
            // std::cout << "==================== " << std::endl;

    const auto &pos = atom_data->atom_struct_owned.position;

    const auto &nlist = neighborlist->neighlist;
    auto &mol_index = atom_data->atom_struct_owned.molecule_index;

    double virialLocal = 0;
    double virialExternalForceLocal = 0;

    int virial_count = 0;
    double r_max = 0;
    double f_max = 0;

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < nlist.size(); ++i)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      const auto type_i = atom_data->atom_struct_owned.type[i];
      const auto mass_inv_i = atom_data->atom_type_params.mass_inv[type_i];
      const auto charge_i = atom_data->atom_type_params.charge[type_i];
      const auto mol_index_i = mol_index[i];
      int id_i = atom_data->atom_struct_owned.id[i];

      for (auto j : nlist[i])
      {
        bool is_ghost = j >= nlist.size();
        Vector<Real_t> pos_j;
        Real_t type_j, mass_inv_j, charge_j;
        int id_j;
        if (is_ghost)
        {
          j -= nlist.size();
          pos_j = atom_data->atom_struct_ghost.position[j];
          type_j = atom_data->atom_struct_ghost.type[j];
          id_j = atom_data->atom_struct_ghost.id[j];
        }
        else
        {
          pos_j = atom_data->atom_struct_owned.position[j];
          type_j = atom_data->atom_struct_owned.type[j];
          id_j = atom_data->atom_struct_owned.id[j];
        }
        charge_j = atom_data->atom_type_params.charge[type_j];
        mass_inv_j = atom_data->atom_type_params.mass_inv[type_j];

        const auto mol_index_j = mol_index[j];

        const auto dr = pos_j - pos[i];

        const auto dr_sq = dr * dr;

        if ((mol_index_i == -1) || (mol_index_i != mol_index_j))
        {
          double forceCoef = 0;
          // virialLocal += - forceCoef * dr_sq;
          // virial_count++;
          // if (r_max < dr_sq) r_max = dr_sq;
          // if (f_max < std::abs(forceCoef)) f_max = std::abs(forceCoef);
          //  if (i == 19)
          //  {
          //    std::cout << "A i=19 , j=" << j<< " dr:" << dr << " drsq: " << dr_sq << " forceCoef: " << forceCoef << " jghost:" <<is_ghost << std::endl;
          //  }
          if (i == 19)
          {
            // std::cout << "XA i=19 , j=" << j<< " dr:" << dr << " drsq: " << dr_sq << " forceCoef: " << forceCoef << " jghost:" <<is_ghost << std::endl;
          }
          if (j == 19)
          {
            // std::cout << "XB i=" << i << " , j=19" << " dr:" << dr.x << " " << dr.y << " " << dr.z << " drsq: " << dr_sq << " forceCoef: " << forceCoef << " jghost:" <<is_ghost << "\n";
            // std::cout << "Xpos_i= " << pos[i].x <<" " << pos[i] << " " << pos[i] <<"\n";
            // std::cout << "Xpos_i= " << pos[i].x <<" " << pos[i] << " " << pos[i] <<"\n";
          }
        }
        else
        {
          double forceCoef = 0;

          if (id_i == 19)
          {

            // std::cout << "SA id_i=19 , id_j=" << id_j<< " dr:" << dr << " drsq: " << dr_sq << " forceCoef: " << forceCoef << " jghost:" <<is_ghost << std::endl;
          }
          if (id_j == 19)
          {
            // std::cout << "SB id_i=" << id_i << " , id_j=19" << " dr:" << dr.x << " " << dr.y << " " << dr.z << " drsq: " << dr_sq << " forceCoef: " << forceCoef << " jghost:" <<is_ghost << "\n";

            // std::cout << "pos_i= " << pos[i] <<"\n";
            // std::cout << "pos_j= " << pos_j <<std::endl;
          }
        }
        if (dr_sq >= cutoff_sq)
          continue;
        const auto dr_norm = std::sqrt(dr_sq);

        auto forceCoef = -k_electrostatic * charge_i * charge_j / (dr_sq * dr_norm) * (1.0 - std::pow(dr_norm / cutoff, beta + 1));
        ;
        const auto force_shifted = forceCoef * dr;

        // const auto force = k_electrostatic * charge_i * charge_j * dr / (dr_sq * dr_norm);
        // const auto force_shifted = force * (1.0 - std::pow(dr_norm / cutoff, beta + 1));
        // if (i < j)
        // virialLocal += - forceCoef * dr_sq;
        //if ((mol_index_i == -1) || (mol_index_i != mol_index_j))
        {
          if (id_i < id_j)
          {
            virialLocal += -forceCoef * dr_sq;
            virial_count++;

            if (id_i == 19 && id_j == 20)
            {
              // std::cout << "SA i: " << i << " j: " << j << " id_i=" << id_i << " , id_j=" << id_j << " dr:" << dr << " drsq: " << dr_sq << " forceCoef: " << forceCoef << " jghost:" << is_ghost << std::endl;
              // std::cout << "SA pi:19: " << pos[i] << std::endl;
              // std::cout << "SA pj:20: " << pos_j << std::endl;
              // std::cout << "ZA i: " << i << " j: " << j << " id_i=" << id_i << " id_j=" << id_j << " drsq: " << dr_sq << " forceCoef: " << forceCoef << " jghost:" << is_ghost << std::endl;

            }
            if (id_i == 20 && id_j == 19)
            {
              // std::cout << "SB i: " << i << " j: " << j << " id_i=" << id_i << " , id_j=" << id_j << " dr:" << dr << " drsq: " << dr_sq << " forceCoef: " << forceCoef << " jghost:" << is_ghost << std::endl;
              // std::cout << "SB pi:19: " << pos[i] << std::endl;
              // std::cout << "SB pj:20: " << pos_j << std::endl;
              // std::cout << "ZB i: " << i << " j: " << j << " id_i=" << id_i << " id_j=" << id_j << " drsq: " << dr_sq << " forceCoef: " << forceCoef << " jghost:" << is_ghost << std::endl;

            }

            if (std::abs(forceCoef) > 498.0 && std::abs(forceCoef) < 499.0)
            {
            }
            if (r_max < dr_sq)
              r_max = dr_sq;
            if (f_max < std::abs(forceCoef))
              f_max = std::abs(forceCoef);
          }
        }

        atom_data->atom_struct_owned.acceleration[i] += force_shifted * mass_inv_i;

        if (!is_ghost)
        {
#ifdef CAVIAR_WITH_OPENMP
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[j].x -= force_shifted.x * mass_inv_j;
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[j].y -= force_shifted.y * mass_inv_j;
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[j].z -= force_shifted.z * mass_inv_j;
#else
          atom_data->atom_struct_owned.acceleration[j] -= force_shifted * mass_inv_j;
#endif
        }
      }

      virialExternalForceLocal += -external_field * charge_i * pos[i];

      const auto force = external_field * charge_i;
      atom_data->atom_struct_owned.acceleration[i] += force * mass_inv_i;
    }

    atom_data->virialForce += virialLocal;
    atom_data->virialExternalForce += virialExternalForceLocal;
    // std::cout << "el virialLocal " << virialLocal << std::endl;
    // std::cout << "el r_max " << r_max << std::endl;
    // std::cout << "el f_max " << f_max << std::endl;
    // std::cout << "el virial_count " << virial_count << std::endl;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
