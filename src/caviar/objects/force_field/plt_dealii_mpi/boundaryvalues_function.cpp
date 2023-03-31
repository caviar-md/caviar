
//========================================================================
//
// Copyright (C) 2019 by deal.II authors and
// Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Some part of this file (Solving Laplace equation) is based on the step-40
// tutorial program of the deal.II library, with extensive modifications
// by Morad Biagooi and Ehsan Nedaaee Oskoee.
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

#ifdef CAVIAR_WITH_DEALII_MPI

#include "caviar/objects/force_field/plt_dealii_mpi.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

namespace caviar {

namespace force_field {

namespace plt_dealii_mpi {

// for ewald in dealii_poisson_mpi:
// if We want to change from Cell List to Verlet List, we have to implement a new
// function for DealII, similar to what it does in 'interpolate_boundary_condition'

double BoundaryValues::potential_of_free_charges (const dealii::Point<3> &p) const
{

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

  const Vector<double> r = {p[0], p[1], p[2]};

  double potential = 0.0;

  for (auto &&f : deal_force -> force_field_custom)
    potential += f->potential (r);

  return potential;
#else
  return 0.0;
#endif
}

double BoundaryValues::value (const Point<3> &p, const unsigned int ) const
{

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  return total_potential - potential_of_free_charges (p);
#else
  return 0.0;
#endif
}

} // plt_dealii
} //force_field

} // namespace caviar
#endif
