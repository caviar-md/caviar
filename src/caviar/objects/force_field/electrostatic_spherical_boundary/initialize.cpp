
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
#include "caviar/utility/macro_constants.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  void Electrostatic_spherical_boundary::initialize()
  {

    auto &pos = atom_data->atom_struct_owned.position;
    auto pos_size = pos.size();

    if (uncharged_particles_optimization)
    {
      int num_of_charged = 0;
      for (unsigned int i = 0; i < pos_size; ++i)
      {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
        auto type_i = atom_data->atom_struct_owned.type[i];
        auto charge_i = atom_data->atom_type_params.charge[type_i];
        if (charge_i != 0.0)
          ++num_of_charged;
      }

      image.position.resize(num_of_charged, Vector<double>{0, 0, 0});
      image.charge.resize(num_of_charged, 0);
    }
    else
    {
      if (image.position.size() < pos_size)
      {
        image.position.resize(pos_size, Vector<double>{0, 0, 0});
        image.charge.resize(pos_size, 0);
      }
    }

    calculate_image_charges();
  }

  //==================================================
  //==================================================
  //==================================================

  void Electrostatic_spherical_boundary::calculate_image_charges()
  {

    const auto &pos = atom_data->atom_struct_owned.position;
    const auto rad_sq = radius * radius;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < pos.size(); ++i)
    {
      #ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      auto type_i = atom_data->atom_struct_owned.type[i];
      double charge_i = atom_data->atom_type_params.charge[type_i];

      if (uncharged_particles_optimization)
        if (charge_i == 0.0)
          continue;

      const auto p_rel = pos[i] - center;
      const auto p_sq = p_rel * p_rel;
      const auto p_sq_inv = 1.0 / p_sq;

      if (p_sq == 0)
        continue;
      if (p_sq > rad_sq)
      {
        error->all(FC_FILE_LINE_FUNC, "particle outside sphere is not implemented yet.");
      }
      else if (p_sq == rad_sq)
      {
        error->all(FC_FILE_LINE_FUNC, "particle is on the shell.");
      }
      else
      {

        image.charge[i] = -charge_i * radius * std::sqrt(p_sq_inv);

        image.position[i] = (p_rel * radius * radius * p_sq_inv) + center;
      }
    }
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
