
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

#ifndef CAVIAR_OBJECTS_WRITER_H
#define CAVIAR_OBJECTS_WRITER_H

#include "caviar/utility/objects_common_headers.h"

namespace caviar {

namespace objects {

/**
 * This class is the base class for all the writers.
 * 
 * 
 */
class Writer : public Pointers {
 public:
  Writer (class CAVIAR *);
  virtual ~Writer ( );
  virtual bool read (class caviar::interpreter::Parser *) = 0;
  virtual void initialize();
  virtual void write();
  virtual void write(int64_t); // current time_step
  virtual void write(double); // current time
  virtual void write(int64_t, double); //time_step and time
  virtual void start_new_files(); //add_time_to_previous
  virtual void start_new_files(std::string &); //add_time_to_previous
  virtual void open_files();
  virtual void close_files();
  virtual void generate();
  bool initialized;
  int64_t last_timestep;
  double last_time;
  double dt;
  int my_mpi_rank, mpi_world_size;
  FC_BASE_OBJECT_COMMON_TOOLS
};

} //objects

} // namespace caviar

#endif
