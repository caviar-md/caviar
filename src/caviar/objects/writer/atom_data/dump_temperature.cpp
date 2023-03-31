
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

  void Atom_data::dump_temperature(int64_t i)
  {
    dump_temperature(i, 0.0);
  }

  void Atom_data::dump_temperature(int64_t i, double t)
  {

    double k_e = atom_data->temperature();

    ofs_temperature << i << " " << t << " " << k_e << std::endl;
  }

} // writer

CAVIAR_NAMESPACE_CLOSE
