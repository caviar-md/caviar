
//========================================================================
//
// Copyright (C) 2019 by deal.II authors and
// Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Some part of this file (Solving Laplace equation) is based on the step-6
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

#ifdef CAVIAR_WITH_DEALII

#include "caviar/objects/force_field/plt_dealii.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/domain.h"
#include "caviar/objects/unique/time_function_3d.h"
//#include "caviar/objects/atom_data.h"

#include <cmath>
#include <iomanip>

namespace caviar {
namespace objects {
namespace force_field {

namespace plt_dealii {

// for ewald in dealii_poisson:
// if We want to change from Cell List to Verlet List, we have to implement a new
// function for DealII, similar to what it does in 'interpolate_boundary_condition'

double BoundaryValues::potential_of_free_charges (const dealii::Point<3> &p) const
{
  const Vector<double> r = {p[0], p[1], p[2]};

  double potential = 0.0;
  if (deal_force->position_offset == nullptr)
    for (auto &&f : deal_force -> force_field_custom)
      potential += f->potential (r);
  else
    for (auto &&f : deal_force -> force_field_custom)
      potential += f->potential (r + deal_force->position_offset->current_value);

  return potential;
}

double BoundaryValues::value (const Point<3> &p, const unsigned int ) const
{
#ifdef CAVIAR_WITH_MPI
  double total_potential_of_free_charges = 0;
  double local_potential_of_free_charges = potential_of_free_charges (p);
//  MPI_Barrier (mpi_communicator); // XXX is it nessecary?

  MPI_Allreduce(&local_potential_of_free_charges,
    &total_potential_of_free_charges,
    1, MPI::DOUBLE, MPI_SUM,  MPI::COMM_WORLD);
    
  return total_potential - total_potential_of_free_charges;
#else
  return total_potential - potential_of_free_charges (p);
#endif     
}

} // plt_dealii
} //force_field
} //objects
} // namespace caviar
#endif
