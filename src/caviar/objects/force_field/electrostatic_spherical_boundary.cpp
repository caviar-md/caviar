
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
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Electrostatic_spherical_boundary::Electrostatic_spherical_boundary(CAVIAR *fptr) : Force_field{fptr}, k_electrostatic{1.0}, external_field{Vector<double>{0, 0, 0}}
  {
    FC_OBJECT_INITIALIZE_INFO
    calculated_once = false;
    radius = 1.0;
    voltage = 0.0;
    center = Vector<double>{0, 0, 0};
    uncharged_particles_optimization = false;
  }

  bool Electrostatic_spherical_boundary::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "uncharged_particles_optimization"))
      {
        uncharged_particles_optimization = true;
      }
      else if (string_cmp(t, "radius"))
      {
        GET_OR_CHOOSE_A_REAL(radius, "", "")
        if (radius < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Force field radius have to non-negative.");
      }
      else if (string_cmp(t, "initialize"))
      {
        initialize();
      }
      else if (string_cmp(t, "voltage"))
      {
        GET_OR_CHOOSE_A_REAL(voltage, "", "")
      }
      else if (string_cmp(t, "external_field"))
      {
        GET_OR_CHOOSE_A_REAL_3D_VECTOR(external_field, "", "");
        // if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Epsilon have to be non-negative.");
      }
      else if (string_cmp(t, "k_electrostatic"))
      {
        GET_OR_CHOOSE_A_REAL(k_electrostatic, "", "")
        if (k_electrostatic < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "k_electrostatic has to be non-negative.");
      }
      else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Electrostatic_spherical_boundary::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
  }

  void Electrostatic_spherical_boundary::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS
    initialize();
    // calculate_image_charges(); called in initialize

    // force-field calculations

    const auto &pos = atom_data->owned.position;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < pos.size(); ++i)
    {

      const auto type_i = atom_data->owned.type[i];
      const auto mass_inv_i = atom_data->owned.mass_inv[type_i];
      const auto charge_i = atom_data->owned.charge[type_i];

      // particle-particle electrostatic force
      /*
      for (unsigned int j=i+1;j<pos.size();++j) {
        const auto type_j = atom_data -> owned.type [j] ;
        const auto mass_inv_j = atom_data -> owned.mass_inv [ type_j ];
        const auto charge_j = atom_data -> owned.charge [ type_j ];
        const auto dr = pos[j] - pos[i];
        const auto dr_sq = dr*dr;
        const auto dr_norm = std::sqrt(dr_sq);
        const auto force = k_electrostatic * charge_i * charge_j * dr / (dr_sq*dr_norm);
        atom_data -> owned.acceleration [i] -= force * mass_inv_i;
        atom_data -> owned.acceleration [j] += force * mass_inv_j;
      }*/

      for (unsigned int j = 0; j < image.position.size(); ++j)
      {
        const auto charge_j = image.charge[j];
        const auto dr = image.position[j] - pos[i];
        const auto dr_sq = dr * dr;
        const auto dr_norm = std::sqrt(dr_sq);

        // in cases that image charge is at infinity or at the origin alongside
        // the particle, this might happen.
        if (dr_norm == 0.0)
          continue;

        const auto force = k_electrostatic * charge_i * charge_j * dr / (dr_sq * dr_norm);
        atom_data->owned.acceleration[i] -= force * mass_inv_i;
      }
    }
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
