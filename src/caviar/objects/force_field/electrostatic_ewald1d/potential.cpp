
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

  double Electrostatic_ewald1d::potential(const Vector<double> &r)
  {
    FC_OBJECT_VERIFY_SETTINGS
    return potential_r(r) + potential_k(r);
  }

  double Electrostatic_ewald1d::potential(const int i)
  {
    FC_OBJECT_VERIFY_SETTINGS
    return potential_r(i) + potential_k(i);
  }

  //======= short range

  // using binlist
  double Electrostatic_ewald1d::potential_r(const Vector<double> &r)
  {
    double potential_value = 0;

    const auto &pos = atom_data->owned.position;
    const auto &binlist = neighborlist->binlist;
    const auto &nb = neighborlist->neigh_bin;
    const auto nb_i = neighborlist->neigh_bin_index(r);
    const int pos_size = pos.size();

    const auto pos_i = r;
    const auto sigma_sq = sigma * sigma;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : potential_value)
#endif
    for (unsigned nb_j = 0; nb_j < nb[nb_i].size(); ++nb_j)
    {
      const auto &nb_ij = nb[nb_i][nb_j];

      for (unsigned i = 0; i < binlist[nb_ij.x][nb_ij.y][nb_ij.z].size(); ++i)
      {

        int j = binlist[nb_ij.x][nb_ij.y][nb_ij.z][i];

        bool is_ghost = j >= pos_size;

        Vector<Real_t> pos_j;
        Real_t type_j;
        if (is_ghost)
        {

          j -= pos_size;
          pos_j = atom_data->ghost.position[j];
          type_j = atom_data->ghost.type[j];
        }
        else
        {
          pos_j = atom_data->owned.position[j];
          type_j = atom_data->owned.type[j];
        }

        const auto charge_j = atom_data->owned.charge[type_j];
        const auto r_ij = pos_i - pos_j;

        if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0)
          continue;

        const auto r_ij_sq = r_ij * r_ij;
        const auto d1 = 1.0 / std::sqrt(r_ij_sq);
        const auto d2 = 1.0 / std::sqrt(r_ij_sq + sigma_sq);
        potential_value += charge_j * (d1 - d2);
      }
    }
    return potential_value * k_electrostatic;
    ;
  }

  //  using neighlist
  double Electrostatic_ewald1d::potential_r(const int i)
  {

    // XXX not checked
    error->all("not implemented. needs fixs for neighlist or maybe impossible.");

    double potential_value = 0;
    const auto &pos = atom_data->owned.position;
    const auto &nlist = neighborlist->neighlist;

    const unsigned pos_size = pos.size();

    const auto &pos_i = atom_data->owned.position[i];
    const auto sigma_sq = sigma * sigma;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : potential_value)
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
        pos_j = atom_data->ghost.position[j];
        type_j = atom_data->ghost.type[j];
      }
      else
      {
        pos_j = atom_data->owned.position[j];
        type_j = atom_data->owned.type[j];
      }

      const auto charge_j = atom_data->owned.charge[type_j];
      const auto r_ij = pos_i - pos_j;

      if (r_ij.x == 0 && r_ij.y == 0 && r_ij.z == 0)
        continue;

      const auto r_ij_sq = r_ij * r_ij;
      const auto d1 = 1.0 / std::sqrt(r_ij_sq);
      const auto d2 = 1.0 / std::sqrt(r_ij_sq + sigma_sq);
      potential_value += coef * charge_j * (d1 - d2);
    }

    return potential_value * k_electrostatic;
  }

  //====== long rang

  double Electrostatic_ewald1d::potential_k(const Vector<double> &r)
  {
    double potential_value = 0;
    const auto &pos = atom_data->owned.position;
    const auto lattice_vec_size = lattice_vec.size();
    const auto sigma_sq = sigma * sigma;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : potential_value)
#endif
    for (unsigned int j = 0; j < pos.size(); ++j)
    {
      const auto type_j = atom_data->owned.type[j];
      const auto charge_j = atom_data->owned.charge[type_j];

      double sum = 0;
      for (unsigned int k = 0; k < lattice_vec_size; ++k)
      {
        const auto dr = r - pos[j] + lattice_vec[k];
        const auto dr_sq = dr * dr;
        // if (dr_sq == 0.0) continue; //XXX removing self potential?
        sum += 1.0 / std::sqrt(dr_sq + sigma_sq);
      }

      potential_value += charge_j * sum;
    }
    return potential_value * k_electrostatic;
  }

  double Electrostatic_ewald1d::potential_k(const int i)
  {
    double potential_value = 0;
    const auto &pos = atom_data->owned.position;
    const auto lattice_vec_size = lattice_vec.size();
    const auto sigma_sq = sigma * sigma;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : potential_value)
#endif
    for (unsigned int j = 0; j < pos.size(); ++j)
    {
      const auto type_j = atom_data->owned.type[j];
      const auto charge_j = atom_data->owned.charge[type_j];

      double sum = 0;
      for (unsigned int k = 0; k < lattice_vec_size; ++k)
      {
        const auto dr = pos[i] - pos[j] + lattice_vec[k];
        const auto dr_sq = dr * dr;
        // if (dr_sq == 0.0) continue; //XXX removing self potential?
        sum += 1.0 / std::sqrt(dr_sq + sigma_sq);
      }

      potential_value += charge_j * sum;
    }
    return potential_value * k_electrostatic;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
