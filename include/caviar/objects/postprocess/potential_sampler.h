
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

#ifndef CAVIAR_OBJECTS_POSTPROCESS_POTENTIALSAMPLER_H
#define CAVIAR_OBJECTS_POSTPROCESS_POTENTIALSAMPLER_H

#include "caviar/objects/unique.h"
#include "caviar/utility/vector.h"
#include <vector>

namespace caviar {

namespace objects {
class Atom_data;
class Md_simulator;

namespace unique {
class Grid_1D;
}

namespace postprocess {

/**
 * This class reads an xyz file, re-runs some frames, then sample
 * the potential values in the given 3D grid, and output them.
 */

class Potential_sampler  : public Postprocess {
 public:
  Potential_sampler (class CAVIAR *);    
  
  ~Potential_sampler () ;

  bool read (caviar::interpreter::Parser *);

  std::vector<caviar::Vector<double>> sampling_position;
  std::vector<caviar::Vector<int>> sampling_position_index;

  void run();

  // this function reads the next frame of the input file and 
  // set the coordinates of the atoms into the atom_data
  void read_next_frame(bool set_the_frame);
  
  void sample_potential();
  
  int step_current;
  
  int step_start;
  int step_end;
  int step_increment;
  
  std::string input_xyz_file_name;
  std::string output_file_name;

  class caviar::objects::Atom_data *atom_data;
  class caviar::objects::Md_simulator *md_simulator;
  
  std::vector<objects::Force_field *> force_field; 
  
  class caviar::objects::unique::Grid_1d *grid_x;
  class caviar::objects::unique::Grid_1d *grid_y;
  class caviar::objects::unique::Grid_1d *grid_z;

  std::ofstream ofs_out;
};

} //postproces
} //objects

} // namespace caviar

#endif
