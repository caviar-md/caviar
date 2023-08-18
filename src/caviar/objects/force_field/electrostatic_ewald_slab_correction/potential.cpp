
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

#include "caviar/objects/force_field/electrostatic_ewald_slab_correction.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/macro_constants.h"

#include <complex>
#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  double Electrostatic_ewald_slab_correction::potential(const Vector<double> &r)
  {
    double sum_p = 0;

    // XXX working scheme of order N. Makes the energy of the order N^2;
    //
    const auto &pos = atom_data->atom_struct_owned.position;
    const auto &type = atom_data->atom_struct_owned.type;
    const auto &charge = atom_data->atom_type_params.charge;
    const auto pos_size = pos.size();

    // XXX no OpenMP parallel yet (due to boolean flag)
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ : sum_p)
#endif
    for (unsigned int j = 0; j < pos_size; ++j)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[j] != my_mpi_rank)
        continue;
#endif
      const auto p = give_slab_local_coordinates(r - pos[j]);
      const double q_j = charge[type[j]];

      bool do_iy_loop_once = true;
      double sum_kx = 0;
      double sum_ky = 0;
      double sum_kp = 0;
      int ip = 0;
      for (auto ix = 0; ix < kx_max; ++ix)
      {
        const auto coshkxz = std::cosh(kx[ix] * p.z);
        const auto coskxx = std::cos(kx[ix] * p.x);

        sum_kx += kx_coef[ix] * coshkxz * coskxx;

        for (auto iy = 0; iy < ky_max; ++iy)
        {
          const auto coshkpz = std::cosh(kp[ip] * p.z);
          const auto coskyy = std::cos(ky[iy] * p.y);
          sum_kp += kp_coef[ip] * coshkpz * coskxx * coskyy;
          if (do_iy_loop_once)
          {
            sum_ky += ky_coef[iy] * std::cosh(ky[iy] * p.z) * coskyy;
          }
          ++ip;
        }
        do_iy_loop_once = false;
      }

      sum_p += q_j * ((2.0 * sum_kp) + sum_kx + sum_ky);
    }

    return sum_p * slab_sum_e_coef * k_electrostatic; // XXX add dipole potential here
  }

  double Electrostatic_ewald_slab_correction::potential(const int i)
  {
    return potential(atom_data->atom_struct_owned.position[i]);
  }

  double Electrostatic_ewald_slab_correction::dipole_potential(const Vector<double> &r)
  {
    return 0 * r.x; // XXX
  }

  double Electrostatic_ewald_slab_correction::dipole_potential(const int i)
  {
    return 0 * i; // XXX is it correct? Since the dipole field is position independent
                  // We also expect that the potential be periodic. Do we?
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
