
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
#include "caviar/objects/atom_data.h"

#include <cmath>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Vector<double> Electrostatic_short_range::field(const Vector<double> &r)
  {
    if (!initialized)
      initialize();
    Vector<double> field_shifted_sum{0, 0, 0};
    const auto &pos = atom_data->owned.position;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : field_shifted_sum)
#endif
    for (unsigned int j = 0; j < pos.size(); ++j)
    {

      const auto type_j = atom_data->owned.type[j];
      const auto charge_j = atom_data->owned.charge[type_j];
      const auto dr = r - pos[j];
      const auto dr_sq = dr * dr;
      if (dr_sq == 0.0)
        continue;
      const auto dr_norm = std::sqrt(dr_sq);
      const auto field = charge_j * dr / (dr_sq * dr_norm);
      field_shifted_sum += field * (1.0 - std::pow(dr_norm / cutoff, beta + 1));
    }
    return field_shifted_sum * k_electrostatic;
  }

  Vector<double> Electrostatic_short_range::field(const int i)
  {
    if (!initialized)
      initialize();
    Vector<double> field_shifted_sum{0, 0, 0};
    const auto &pos = atom_data->owned.position;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : field_shifted_sum)
#endif
    for (unsigned int j = 0; j < pos.size(); ++j)
    {
      if (i == static_cast<int>(j))
        continue;
      const auto type_j = atom_data->owned.type[j];
      const auto charge_j = atom_data->owned.charge[type_j];
      const auto dr = pos[i] - pos[j];
      const auto dr_sq = dr * dr;
      const auto dr_norm = std::sqrt(dr_sq);
      const auto field = charge_j * dr / (dr_sq * dr_norm);

      field_shifted_sum += field * (1.0 - std::pow(dr_norm / cutoff, beta + 1));
    }
    return field_shifted_sum * k_electrostatic;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
