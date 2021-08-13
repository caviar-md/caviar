
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

#ifndef CAVIAR_OBJECTS_WRITER_FORCEFIELD_H
#define CAVIAR_OBJECTS_WRITER_FORCEFIELD_H

#include "caviar/objects/writer.h"

namespace caviar {
namespace objects {
namespace writer {

/**
 * This class has a writer for force-field
 * 
 * 
 */
class Force_field : public Writer {
 public:
  Force_field (class CAVIAR *);
  ~Force_field ( );
  bool read (class caviar::interpreter::Parser *);
  void initialize();
  void write();
  void write(int64_t); // current time_step
  void write(double); // current time
  void write(int64_t, double); //time_step and time
  void start_new_files(); //add_time_to_previous
  void start_new_files(std::string &); //add_time_to_previous

 public:

};

} //writer
} //objects
} // namespace caviar

#endif
