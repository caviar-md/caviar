
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

#include "caviar/objects/force_field/electrostatic_ewald_k.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/macro_constants.h"

#include <cmath>
#include <iomanip>

#ifdef CAVIAR_WITH_OPENMP

#pragma omp declare reduction(+                       \
                              : std::complex <double> \
                              : omp_out += omp_in)
#endif

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Vector<double> Electrostatic_ewald_k::k_space_field(const int i)
  {
    Vector<double> field{0, 0, 0};
    const auto &pos = atom_data->atom_struct_owned.position;
    static std::complex<double> ii(0.0, 1.0);

    std::complex<double> sum_kx(0.0, 0.0);
    std::complex<double> sum_ky(0.0, 0.0);
    std::complex<double> sum_kz(0.0, 0.0);

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : sum_kx, sum_ky, sum_kz)
#endif
    for (int k = 0; k < n_k_vectors; ++k)
    {

      const auto sum_j_c = std::conj(potential_k_coef_cmplx[k]);
      const auto c = field_k_coef[k] * sum_j_c * std::exp(ii * (k_vector[k] * pos[i]));

      sum_kx += k_vector[k].x * c;
      sum_ky += k_vector[k].y * c;
      sum_kz += k_vector[k].z * c;
    }

    const double sum_jx = std::imag(sum_kx);
    const double sum_jy = std::imag(sum_ky);
    const double sum_jz = std::imag(sum_kz);

    Vector<double> sv{sum_jx, sum_jy, sum_jz};

    const auto sum = FC_4PI * sv;

    field = k_electrostatic * l_xyz_inv * sum;

    return field;
  }

  Vector<double> Electrostatic_ewald_k::k_space_field(const Vector<double> &r)
  {
    Vector<double> field{0, 0, 0};

    static std::complex<double> ii(0.0, 1.0);

    std::complex<double> sum_kx(0.0, 0.0);
    std::complex<double> sum_ky(0.0, 0.0);
    std::complex<double> sum_kz(0.0, 0.0);

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : sum_kx, sum_ky, sum_kz)
#endif
    for (int k = 0; k < n_k_vectors; ++k)
    {

      const auto sum_j_c = std::conj(potential_k_coef_cmplx[k]);
      const auto c = field_k_coef[k] * sum_j_c * std::exp(ii * (k_vector[k] * r));

      sum_kx += k_vector[k].x * c;
      sum_ky += k_vector[k].y * c;
      sum_kz += k_vector[k].z * c;
    }

    const double sum_jx = std::imag(sum_kx);
    const double sum_jy = std::imag(sum_ky);
    const double sum_jz = std::imag(sum_kz);

    Vector<double> sv{sum_jx, sum_jy, sum_jz};

    const auto sum = FC_4PI * sv;

    field = k_electrostatic * l_xyz_inv * sum;

    return field;
  }

  Vector<double> Electrostatic_ewald_k::dipole_field()
  {
    return dipole_field_vector;
  }

  Vector<double> Electrostatic_ewald_k::field(int i)
  {
    auto field = k_space_field(i);
    if (dipole)
      field += dipole_field_vector;
    return field;
  }

  Vector<double> Electrostatic_ewald_k::field(const Vector<double> &r)
  {
    auto field = k_space_field(r);
    if (dipole)
      field += dipole_field_vector;
    return field;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
