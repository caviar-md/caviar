
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

#include <cmath>
#include <iomanip>
#include <complex>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Electrostatic_ewald_slab_correction::Electrostatic_ewald_slab_correction(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    slab_normal_axis = 2;

    initialized = false;
    calculated_once = false;
    k_electrostatic = 1.0;
    kx_max = 1;
    ky_max = 1;
    lz = 1.0;
  }

  bool Electrostatic_ewald_slab_correction::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "k_max"))
      {
        int k_max = 1;
        GET_OR_CHOOSE_A_INT(k_max, "", "")
        if (k_max < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "'k_max' have to non-negative.");
        kx_max = ky_max = k_max;
      }
      else if (string_cmp(t, "kx_max"))
      {
        GET_OR_CHOOSE_A_INT(kx_max, "", "")
        if (kx_max < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "'kx_max' have to non-negative.");
      }
      else if (string_cmp(t, "ky_max"))
      {
        GET_OR_CHOOSE_A_INT(ky_max, "", "")
        if (ky_max < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "'ky_max' have to non-negative.");
      }
      else if (string_cmp(t, "lz"))
      {
        GET_OR_CHOOSE_A_REAL(lz, "", "")
        if (lz < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "'lz' has to be non-negative.");
      }
      else if (string_cmp(t, "slab_normal_axis"))
      {
        const auto t = parser->get_val_token();
        if (t.kind == caviar::interpreter::Kind::identifier)
        {
          const auto ts = t.string_value;
          if (string_cmp(ts, "x"))
            slab_normal_axis = 0;
          if (string_cmp(ts, "y"))
            slab_normal_axis = 1;
          if (string_cmp(ts, "z"))
            slab_normal_axis = 2;
        }
        else
          error->all(FC_FILE_LINE_FUNC_PARSE, "expected 'x', 'y' or 'z'.");
        if (slab_normal_axis == 0)
          output->info("slab local axis 'x' = simulation "
                       "global 'y' axis, slab local axis 'y' = simulation global 'z' axis.");

        if (slab_normal_axis == 1)
          output->info("slab local axis 'x' = simulation "
                       "global 'x' axis, slab local axis 'y' = simulation global 'z' axis.");

        if (slab_normal_axis == 2)
          output->info("slab local axis 'x' = simulation "
                       "global 'x' axis, slab local axis 'y' = simulation global 'y' axis.");
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
      else if (string_cmp(t, "set_domain") || string_cmp(t, "domain"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Electrostatic_ewald_slab_correction::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(domain)
    my_mpi_rank = atom_data->get_mpi_rank();
  }

  void Electrostatic_ewald_slab_correction::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    initialize();

    // XXX Scheme using field functions. Of the order of the field implementation.
    // /*
    const auto &pos = atom_data->atom_struct_owned.position;
    const unsigned pos_size = pos.size();
    double virialLocal = 0;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < pos_size; ++i)
    {
      #ifdef CAVIAR_WITH_MPI
    if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;
#endif

      const auto pos_i = atom_data->atom_struct_owned.position[i];
      const auto type_i = atom_data->atom_struct_owned.type[i];
      const auto charge_i = atom_data->atom_type_params.charge[type_i];
      const auto mass_inv_i = atom_data->atom_type_params.mass_inv[type_i];

      const auto force = charge_i * field(pos_i); //
                                                  //    const auto force = charge_i * field (i);  //
      atom_data->atom_struct_owned.acceleration[i] += force * mass_inv_i;
    }

    atom_data->virialForce += virialLocal;

    // */
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
