
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

#include "caviar/objects/writer.h"
#include "caviar/interpreter/communicator.h"

namespace caviar {

namespace objects {

Writer::Writer (CAVIAR *fptr) : Pointers{fptr}, initialized{false}, my_mpi_rank{comm->me},
    mpi_world_size{comm->nprocs}
 {
  FC_OBJECT_INITIALIZE
}

Writer::~Writer () {}

void Writer::verify_settings () {
  
}

void Writer::initialize(){}
void Writer::write(){}
void Writer::write(int64_t){} // current time_step
void Writer::write(double){} // current time
void Writer::write(int64_t, double){} //time_step and time
void Writer::start_new_files(){} //add_time_to_previous
void Writer::start_new_files(std::string &){} //add_time_to_previous
void Writer::open_files() {}
void Writer::close_files() {}
void Writer::generate() {}
} //objects


} // namespace caviar

