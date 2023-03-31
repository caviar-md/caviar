
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

#include "caviar/objects/postprocess.h"
#include "caviar/utility/vector.h"
#include <vector>

namespace caviar {


class Atom_data;
class Md_simulator;
class Force_field;

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
  int read_next_frame(bool set_the_frame, bool read_velocity);
  
  void sample_potential();
  
  void set_positions_vectors();
  
  int step_current;
  
  int step_start;
  int step_end;
  int step_increment;
  
  bool read_velocity;
  
  std::string input_xyz_file_name;
  std::string output_file_name;

  class caviar::Atom_data *atom_data;
  class caviar::Md_simulator *md_simulator;
  
  std::vector<Force_field *> force_field; 
  
  class caviar::unique::Grid_1D *grid_x;
  class caviar::unique::Grid_1D *grid_y;
  class caviar::unique::Grid_1D *grid_z;

  std::ofstream ofs_out;
};

} //postproces


} // namespace caviar

#endif
