
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

#include "caviar/objects/writer/force_field.h"
#include "caviar/utility/interpreter_io_headers.h"
// #include <ctime>

CAVIAR_NAMESPACE_OPEN

namespace writer
{
  /*
  void Md_simulator::output_total_force() {
    FC_NULLPTR_CHECK(atom_data)
    std::cout << std::setprecision(10);
    std::cout << "Md_simulator::output_total_force :\n";
    const auto &pos_size = atom_data->owned.position.size();
    const auto &acc = atom_data->owned.acceleration;
    const auto &type = atom_data->owned.type;
    const auto &mass = atom_data->owned.mass;

    for (unsigned i = 0; i < pos_size ; ++i)
      std::cout << i << ": " << acc[i]*mass[ type[i] ] << "\n";
    std::cout << std::setprecision(6);
  }

  void Md_simulator::output_total_energy() {
    if (force_field.size()==0) {
      std::cout << "Md_simulator::output_total_energy: force_field.size()=0 \n";
      return;
    }
    std::cout << std::setprecision(10);
    double sum_e = 0;
    for (auto f : force_field) {
      const auto e = f -> energy ();
      sum_e += e;
      std::cout << "e: " << e << "\n";
    }
    std::cout << "sum_e: " << sum_e << std::endl;
    std::cout << std::setprecision(6);
  }
  */
} // Force_field

CAVIAR_NAMESPACE_CLOSE
