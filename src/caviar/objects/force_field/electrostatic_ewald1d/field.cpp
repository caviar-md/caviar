
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

#include "caviar/objects/force_field/electrostatic_ewald1d.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/neighborlist.h"

#include <cmath>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  //======= total potential

  Vector<double> Electrostatic_ewald1d::field(const Vector<double> &r)
  {
    FC_OBJECT_VERIFY_SETTINGS
    return field_r(r) + field_k(r);
  }

  Vector<double> Electrostatic_ewald1d::field(const int i)
  {
    FC_OBJECT_VERIFY_SETTINGS
    return field_r(i) + field_k(i);
  }

  //======= short range

  Vector<double> Electrostatic_ewald1d::field_r(const Vector<double> &r)
  {

    Vector<double> field{0, 0, 0};

    const auto &pos = atom_data->atom_struct_owned.position;
    const unsigned pos_size = pos.size();

    const auto pos_i = r;

    const auto &binlist = neighborlist->binlist;
    const auto &nb = neighborlist->neigh_bin;
    const auto nb_i = neighborlist->neigh_bin_index(r);
    const auto sigma_sq = sigma * sigma;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : field)
#endif
    for (unsigned nb_j = 0; nb_j < nb[nb_i].size(); ++nb_j)
    {
      const auto &nb_ij = nb[nb_i][nb_j];

      for (unsigned i = 0; i < binlist[nb_ij.x][nb_ij.y][nb_ij.z].size(); ++i)
      {

        unsigned int j = binlist[nb_ij.x][nb_ij.y][nb_ij.z][i];

        bool is_ghost = j >= pos_size;
        Vector<Real_t> pos_j;
        Real_t type_j;
        if (is_ghost)
        {
          j -= pos_size;
          pos_j = atom_data->atom_struct_ghost.position[j];
          type_j = atom_data->atom_struct_ghost.type[j];
        }
        else
        {
          pos_j = atom_data->atom_struct_owned.position[j];
          type_j = atom_data->atom_struct_owned.type[j];
        }

        const auto charge_j = atom_data->atom_type_params.charge[type_j];

        const auto r_ij = pos_i - pos_j;

        if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0)
          continue;

        const auto r_ij_sq = r_ij * r_ij;
        const auto d1 = 1.0 / std::sqrt(r_ij_sq);
        const auto d2 = 1.0 / std::sqrt(r_ij_sq + sigma_sq);

        const auto sum_r = r_ij * ((d1 * d1 * d1) - (d2 * d2 * d2));

        field += charge_j * sum_r;
      }
    }

    return field * k_electrostatic;
    ;
  }

  Vector<double> Electrostatic_ewald1d::field_r(const int i)
  {
    Vector<double> field{0, 0, 0};

    error->all("not implemented. needs fixs for neighlist or maybe impossible.");

    const auto &pos = atom_data->atom_struct_owned.position;
    const unsigned pos_size = pos.size();

    const auto &nlist = neighborlist->neighlist;

    const auto pos_i = atom_data->atom_struct_owned.position[i];
    const auto sigma_sq = sigma * sigma;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : field)
#endif
    // for (auto j : nlist[i]) {
    for (unsigned int k = 0; k < nlist[i].size(); ++k)
    {
      auto j = nlist[i][k];
      double coef = 2.0; // ewald: 'coef=2' for owned in 'neighlist'. Not for binlist.
      bool is_ghost = j >= pos_size;
      Vector<Real_t> pos_j;
      Real_t type_j;
      if (is_ghost)
      {
        coef = 1.0; // ewald:'coef=1' for ghost in 'neighlist'. Not for binlist.
        j -= pos_size;
        pos_j = atom_data->atom_struct_ghost.position[j];
        type_j = atom_data->atom_struct_ghost.type[j];
      }
      else
      {
        pos_j = atom_data->atom_struct_owned.position[j];
        type_j = atom_data->atom_struct_owned.type[j];
      }

      const auto charge_j = atom_data->atom_type_params.charge[type_j];

      const auto r_ij = pos_i - pos_j;

      if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0)
        continue;

      const auto r_ij_sq = r_ij * r_ij;
      const auto d1 = 1.0 / std::sqrt(r_ij_sq);
      const auto d2 = 1.0 / std::sqrt(r_ij_sq + sigma_sq);
      const auto sum_r = r_ij * ((d1 * d1 * d1) - (d2 * d2 * d2));

      field += coef * charge_j * sum_r;
    }

    return field * k_electrostatic;
  }

  //======= long range

  Vector<double> Electrostatic_ewald1d::field_k(const Vector<double> &r)
  {
    Vector<double> field{0, 0, 0};

    const auto &pos = atom_data->atom_struct_owned.position;
    const auto lattice_vec_size = lattice_vec.size();
    const auto sigma_sq = sigma * sigma;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : field)
#endif
    for (unsigned int j = 0; j < pos.size(); ++j)
    {
      const auto type_j = atom_data->atom_struct_owned.type[j];
      const auto charge_j = atom_data->atom_type_params.charge[type_j];

      Vector<double> sum{0, 0, 0};
      for (unsigned int k = 0; k < lattice_vec_size; ++k)
      {
        const auto dr = r - pos[j] + lattice_vec[k];
        const auto dr_sq = dr * dr;
        // if (dr_sq == 0.0) continue; //XXX removing self field?
        const auto d2 = 1.0 / std::sqrt(dr_sq + sigma_sq);
        sum += dr * (d2 * d2 * d2);
      }

      field += charge_j * sum;
    }

    return field * k_electrostatic;
  }

  Vector<double> Electrostatic_ewald1d::field_k(const int i)
  {
    Vector<double> field{0, 0, 0};

    const auto &pos = atom_data->atom_struct_owned.position;
    const auto lattice_vec_size = lattice_vec.size();
    const auto sigma_sq = sigma * sigma;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : field)
#endif
    for (unsigned int j = 0; j < pos.size(); ++j)
    {
      const auto type_j = atom_data->atom_struct_owned.type[j];
      const auto charge_j = atom_data->atom_type_params.charge[type_j];

      Vector<double> sum{0, 0, 0};
      for (unsigned int k = 0; k < lattice_vec_size; ++k)
      {
        const auto dr = pos[i] - pos[j] + lattice_vec[k];
        const auto dr_sq = dr * dr;
        // if (dr_sq == 0.0) continue; //XXX removing self field?
        const auto d2 = 1.0 / std::sqrt(dr_sq + sigma_sq);
        sum += dr * (d2 * d2 * d2);
      }

      field += charge_j * sum;
    }

    return field * k_electrostatic;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
