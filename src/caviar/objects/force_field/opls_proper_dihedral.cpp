
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Implementation of Opls_proper_dihedral is done by Nasrin Eyvazi
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

#include "caviar/objects/force_field/opls_proper_dihedral.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

/// U(teta) = 1/2 {K1[1+cos(teta)] + K2[1-cos(2*teta)] + K3[1+cos(3*teta)] - K4[1-cos(4*teta)] }
CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Opls_proper_dihedral::Opls_proper_dihedral(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
  }

  bool Opls_proper_dihedral::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "dihedral_coef1"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(dihedral_coef1)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Dihedral coef1. have to be non-negative.");
        //}  else if (string_cmp(t,"dissip_coef")) {
        //  GET_A_STDVECTOR_REAL_ELEMENT(dissip_coef)
        //  if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Dissipation coef. have to be non-negative.");
      }
      else if (string_cmp(t, "dihedral_coef2"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(dihedral_coef2)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Dihedral coef2. have to be non-negative.");
      }
      else if (string_cmp(t, "dihedral_coef3"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(dihedral_coef3)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Dihedral coef3. have to be non-negative.");
      }
      else if (string_cmp(t, "dihedral_coef4"))
      {
        GET_A_STDVECTOR_REAL_ELEMENT(dihedral_coef4)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Dihedral coef4. have to be non-negative.");
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

  void Opls_proper_dihedral::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(domain)
    my_mpi_rank = atom_data->get_mpi_rank();
  }

  void Opls_proper_dihedral::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    auto &pos = atom_data->atom_struct_owned.position;
    auto &type = atom_data->atom_struct_owned.type;
    auto &mass_inv = atom_data->atom_type_params.mass_inv;


#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < atom_data->molecule_struct_owned.size(); i++)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->molecule_struct_owned[i].ghost)
        continue;
#endif
    auto &atomic_properdihedral_vector = atom_data->molecule_struct_owned[i].atomic_properdihedral_vector;

      for (unsigned int j = 0; j < atomic_properdihedral_vector.size(); j++)
      {
        int id_1 = atomic_properdihedral_vector[j].id_1;
        int id_2 = atomic_properdihedral_vector[j].id_2;
        int id_3 = atomic_properdihedral_vector[j].id_3;
        int id_4 = atomic_properdihedral_vector[j].id_4;

        int k1 = atom_data->atom_id_to_index[id_1];
        int k2 = atom_data->atom_id_to_index[id_2];
        int k3 = atom_data->atom_id_to_index[id_3]; 
        int k4 = atom_data->atom_id_to_index[id_4];

        int dtype = atomic_properdihedral_vector[j].type;

#if defined(CAVIAR_WITH_MPI)
        auto p21 = pos[k1] - pos[k2];
        auto p23 = pos[k3] - pos[k2];
        auto o3 = pos[k3] - ((pos[k2] + pos[k3]) / 2);
        auto p34 = pos[k4] - pos[k3];
#else
        auto p21 = domain->periodic_distance(pos[k1] - pos[k2]);
        auto p23 = domain->periodic_distance(pos[k3] - pos[k2]);
        auto o3 = domain->periodic_distance(pos[k3] - ((pos[k2] + pos[k3]) / 2));
        auto p34 = domain->periodic_distance(pos[k4] - pos[k3]);
#endif

        auto m = cross_product(p21, p23);
        auto n = cross_product(p23, p34);
        auto m_size_inv = 1.0 / norm(m);
        auto n_size_inv = 1.0 / norm(n);

        auto teta_cos = (m * n) * (m_size_inv * n_size_inv);
        auto teta = std::acos(teta_cos);

        auto torque = 0.5 * (dihedral_coef1[dtype] * sin(teta) - 2 * dihedral_coef2[dtype] * sin(2 * teta) + 3 * dihedral_coef3[dtype] * sin(3 * teta) - 4 * dihedral_coef4[dtype] * sin(4 * teta));

        auto p21_size_inv = 1.0 / norm(p21);
        auto p23_size_inv = 1.0 / norm(p23);
        auto p34_size_inv = 1.0 / norm(p34);
        auto o3_size_inv = 1.0 / norm(o3);
        auto p21_norm = p21 * p21_size_inv;
        auto p23_norm = p23 * p23_size_inv;
        auto p34_norm = p34 * p34_size_inv;

        auto f1 = cross_product(p21_norm, p23_norm);
        auto f2 = cross_product(p23_norm, p34_norm);
        auto angle123_cos = (p21 * p23) * (p21_size_inv * p23_size_inv);
        auto angle123 = std::acos(angle123_cos);
        auto angle234_cos = ((-p23) * p34) * (p23_size_inv * p34_size_inv);
        auto angle234 = std::acos(angle234_cos);

        auto force_1 = torque * f1 * (p21_size_inv / sin(angle123));
        auto force_4 = torque * f2 * (p34_size_inv / sin(angle234));

        auto tc = -(cross_product(o3, force_4) + 0.5 * cross_product(p34, force_4) + 0.5 * cross_product(p21, force_1));

        auto force_3 = (o3_size_inv) * (o3_size_inv)*cross_product(tc, o3);

        {
#ifdef CAVIAR_WITH_OPENMP
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k1].x += force_1.x * mass_inv[type[k1]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k1].y += force_1.y * mass_inv[type[k1]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k1].z += force_1.z * mass_inv[type[k1]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k4].x += force_4.x * mass_inv[type[k4]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k4].y += force_4.y * mass_inv[type[k4]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k4].z += force_4.z * mass_inv[type[k4]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k3].x += force_3.x * mass_inv[type[k3]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k3].y += force_3.y * mass_inv[type[k3]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k3].z += force_3.z * mass_inv[type[k3]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k2].x -= (force_1.x + force_3.x + force_4.x) * mass_inv[type[k2]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k2].y -= (force_1.y + force_3.y + force_4.y) * mass_inv[type[k2]];
#pragma omp atomic
          atom_data->atom_struct_owned.acceleration[k2].z -= (force_1.z + force_3.z + force_4.z) * mass_inv[type[k2]];
#else
          atom_data->atom_struct_owned.acceleration[k1] += force_1 * mass_inv[type[k1]];
          atom_data->atom_struct_owned.acceleration[k4] += force_4 * mass_inv[type[k4]];
          atom_data->atom_struct_owned.acceleration[k3] += force_3 * mass_inv[type[k3]];
          atom_data->atom_struct_owned.acceleration[k2] -= (force_1 + force_3 + force_4) * mass_inv[type[k2]];
#endif
        }
      }
    }
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
