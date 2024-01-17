
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

#include "caviar/objects/force_field/spring_angle.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Spring_angle::Spring_angle(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
  }

  bool Spring_angle::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR
      if (string_cmp(t, "elastic_coef"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(elastic_coef)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Elastic coef. have to be non-negative.");
      }
      else if (string_cmp(t, "dissip_coef"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(dissip_coef)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Dissipation coef. have to be non-negative.");
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

  void Spring_angle::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(domain)
    my_mpi_rank = atom_data->get_mpi_rank();

    unsigned int btypeMax = 0;
    for (unsigned int i = 0; i < atom_data->molecule_struct_owned.size(); i++)
    {
      auto &atomic_angle_vector = atom_data->molecule_struct_owned[i].atomic_angle_vector;

      for (unsigned int j = 0; j < atomic_angle_vector.size(); j++)
      {
        unsigned int btype = atomic_angle_vector[j].type;
        if (btypeMax < btype)
          btypeMax = btype;
      }
    }

    if (elastic_coef.size() < btypeMax + 1)
    {
      elastic_coef.resize(btypeMax + 1, 0);
      output->warning("Spring_angle::verify_settings: resizing elastic_coef");
    }

    if (dissip_coef.size() < btypeMax + 1)
    {
      dissip_coef.resize(btypeMax + 1, 0);
      output->warning("Spring_angle::verify_settings: resizing dissip_coef");
    }
  }

  void Spring_angle::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    auto &pos = atom_data->atom_struct_owned.position;
    // auto &pos_old = atom_data -> atom_struct_owned.position_old;
    // auto &vel = atom_data -> atom_struct_owned.velocity;
    auto &type = atom_data->atom_struct_owned.type;
    auto &mass_inv = atom_data->atom_type_params.mass_inv;
    double virialLocal = 0;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < atom_data->molecule_struct_owned.size(); i++)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->molecule_struct_owned[i].ghost)
        continue;
#endif

      auto &atomic_angle_vector = atom_data->molecule_struct_owned[i].atomic_angle_vector;

      for (unsigned int j = 0; j < atomic_angle_vector.size(); j++)
      {
        int atype = atomic_angle_vector[j].type;

        if (elastic_coef[atype] == 0 && dissip_coef[atype] == 0)
          continue;

        int id_1 = atomic_angle_vector[j].id_1;
        int id_2 = atomic_angle_vector[j].id_2;
        int id_3 = atomic_angle_vector[j].id_3;

        int k1 = atom_data->atom_id_to_index[id_1], k2 = atom_data->atom_id_to_index[id_2], k3 = atom_data->atom_id_to_index[id_3];



        double angle_value = atomic_angle_vector[j].value;

#if defined(CAVIAR_WITH_MPI)
        auto p21 = pos[k1] - pos[k2];
        auto p23 = pos[k3] - pos[k2];
#else
        auto p21 = domain->periodic_distance(pos[k1] - pos[k2]);
        auto p23 = domain->periodic_distance(pos[k3] - pos[k2]);
#endif
        auto p21_size_inv = 1.0 / norm(p21);
        auto p23_size_inv = 1.0 / norm(p23);

        auto angle_cos = (p21 * p23) * (p21_size_inv * p23_size_inv);
        auto angle = std::acos(angle_cos);
        auto angle_diff = angle - angle_value;

        /*      // this is not a good way.
          double angle_dot = 0.0;
          if (!dissip_coef[atype]==0) {

        #if defined(CAVIAR_WITH_MPI)
            auto p21 = pos_old[k1] - pos_old[k2];
            auto p23 = pos_old[k3] - pos_old[k2];
        #else
            auto p21 = domain-> periodic_distance(pos_old[k1] - pos_old[k2]);
            auto p23 = domain-> periodic_distance(pos_old[k3] - pos_old[k2]);
        #endif
            auto p21_size_inv = 1.0/norm(p21);
            auto p23_size_inv = 1.0/norm(p23);

            auto angle_cos_o = (p21*p23)*(p21_size_inv*p23_size_inv);
            auto angle_o = std::acos(angle_cos);
            auto dt = 0.001;
            angle_dot = (angle - angle_o)/dt;
          }

          auto torque = elastic_coef[atype]*angle_diff - dissip_coef[atype]*angle_dot;
        */

        auto torque = elastic_coef[atype] * angle_diff;

        auto p21_norm = p21 * p21_size_inv;
        auto p23_norm = p23 * p23_size_inv;

        auto n = cross_product(p21_norm, p23_norm);
        auto f12 = cross_product(n, p21_norm);
        auto f32 = cross_product(p23_norm, n);

        // r * F = Tau -> F = Tau / r
        auto force_12 = torque * f12 * p21_size_inv;
        auto force_32 = torque * f32 * p23_size_inv;

#ifdef CAVIAR_WITH_OPENMP
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k1].x += force_12.x * mass_inv[type[k1]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k1].y += force_12.y * mass_inv[type[k1]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k1].z += force_12.z * mass_inv[type[k1]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k3].x += force_32.x * mass_inv[type[k3]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k3].y += force_32.y * mass_inv[type[k3]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k3].z += force_32.z * mass_inv[type[k3]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k2].x -= (force_12.x + force_32.x) * mass_inv[type[k2]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k2].y -= (force_12.y + force_32.y) * mass_inv[type[k2]];
#pragma omp atomic
        atom_data->atom_struct_owned.acceleration[k2].z -= (force_12.z + force_32.z) * mass_inv[type[k2]];
#else
        atom_data->atom_struct_owned.acceleration[k1] += force_12 * mass_inv[type[k1]];
        atom_data->atom_struct_owned.acceleration[k3] += force_32 * mass_inv[type[k3]];
        atom_data->atom_struct_owned.acceleration[k2] -= (force_12 + force_32) * mass_inv[type[k2]];
#endif
      }
    }

    atom_data->virialForce += virialLocal;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
