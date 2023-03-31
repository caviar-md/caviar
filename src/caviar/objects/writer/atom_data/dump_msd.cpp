
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
#include "caviar/objects/domain.h"
#include "caviar/utility/interpreter_io_headers.h"


namespace caviar {



namespace writer {


void Atom_data::dump_msd (int64_t i) {
  dump_msd(i, 0.0);
}

void Atom_data::dump_msd (int64_t i, double t) {

  const auto &pos = atom_data -> owned.position;
  const auto &type = atom_data -> owned.type;

  auto pos_size = pos.size();
  if (i < msd_initial_step) return;
  if (i == msd_initial_step) {
    msd_initial_position.resize(pos_size);
    for (unsigned int j = 0; j < pos_size; ++j) {
      msd_initial_position[j] = pos[j];
      atom_data -> owned.msd_domain_cross[j] = Vector<int> {0,0,0};
    }
    if (my_mpi_rank==0)
      ofs_msd << i << " " << t << " " << "0.0" << "\n";
    return;
  }
  
  if (pos_size != msd_initial_position.size())
    error->all(FC_FILE_LINE_FUNC,"  (pos.size != msd_initial_position.size())");

  // XXX note that this STATIC value may affect NPT ensembles
  static caviar::Vector<double> domain_dx= {(domain->upper_global.x -  domain->lower_global.x),
                                               (domain->upper_global.y -  domain->lower_global.y),
                                               (domain->upper_global.z -  domain->lower_global.z)};

  double sum_dr_sq = 0.0;
  for (unsigned int j = 0; j < pos_size; ++j) {
    if (type[j] != static_cast<unsigned int>(msd_type)) continue;
    auto dr = msd_initial_position[j] - pos[j];

    dr.x += atom_data -> owned.msd_domain_cross[j].x*domain_dx.x;
    dr.y += atom_data -> owned.msd_domain_cross[j].y*domain_dx.y;
    dr.z += atom_data -> owned.msd_domain_cross[j].z*domain_dx.z;

    //domain-> periodic_distance(dr);

    auto dr_sq = dr*dr;
    sum_dr_sq += dr_sq;

  }

  double msd = sum_dr_sq / pos_size;

  if (my_mpi_rank==0) {
    ofs_msd << i << " " << t << " " << msd << std::endl;
  }

}


} // writer 
 

} // namespace caviar


