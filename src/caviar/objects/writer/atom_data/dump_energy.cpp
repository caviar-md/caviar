
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

#include "caviar/objects/writer/atom_data.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace writer
{

  void Atom_data::dump_energy_mpi(int64_t , double )
  {
    // double k_e = atom_data->kinetic_energy();

    // ofs_energy << i << " " << t << " " << k_e << std::endl;
  }

  void Atom_data::dump_energy_mpi_shared_atoms(int64_t i, double t)
  {

    double k_e = atom_data->kinetic_energy();

    if (my_mpi_rank != 0)
      return;

    ofs_energy << i << " " << t << " " << k_e << std::endl;
  }

  void Atom_data::dump_energy_serial(int64_t i, double t, bool)
  {

    double k_e = atom_data->kinetic_energy();

    ofs_energy << i << " " << t << " " << k_e << std::endl;
  }

  void Atom_data::dump_temperature_mpi(int64_t , double )
  {
    // double k_e = atom_data->temperature();

    // ofs_temperature << i << " " << t << " " << k_e << std::endl;
  }

  void Atom_data::dump_temperature_mpi_shared_atoms(int64_t i, double t)
  {
    double k_e = atom_data->temperature();

    if (my_mpi_rank != 0)
      return;

    ofs_temperature << i << " " << t << " " << k_e << std::endl;
  }

  void Atom_data::dump_temperature_serial(int64_t i, double t, bool )
  {

    double k_e = atom_data->temperature();

    ofs_temperature << i << " " << t << " " << k_e << std::endl;
  }

  void Atom_data::dump_pressure_mpi(int64_t , double )
  {
    // double p = atom_data->pressure();

    // ofs_pressure << i << " " << t << " " << p << std::endl;
  }

  void Atom_data::dump_pressure_mpi_shared_atoms(int64_t i, double t)
  {

    double p = atom_data->pressure();

    if (my_mpi_rank != 0)
      return;

    ofs_pressure << i << " " << t << " " << p << std::endl;
  }

  void Atom_data::dump_pressure_serial(int64_t i, double t, bool )
  {

    double p = atom_data->pressure();

    ofs_pressure << i << " " << t << " " << p << std::endl;
  }

} // writer

CAVIAR_NAMESPACE_CLOSE
