
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

#include "caviar/objects/force_field/electrostatic_spherical_boundary.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"

#include <cmath>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  double Electrostatic_spherical_boundary::energy()
  {
    // XXX scheme using potential formula.
    const auto &pos = atom_data->atom_struct_owned.position;
    double energy_r = 0;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ \
                                   : energy_r)
#endif
    for (unsigned int j = 0; j < pos.size(); ++j)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[j] != my_mpi_rank)
        continue;
#endif
      const auto type_j = atom_data->atom_struct_owned.type[j];
      const auto charge_j = atom_data->atom_type_params.charge[type_j];
      // energy_r += charge_j * potential(j); //
      energy_r += charge_j * potential(pos[j]); //
    }

    // there's no particle-particle interaction. only image-particle. so there
    // will be no '0.5' coefficient for energy.
    return energy_r;
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
