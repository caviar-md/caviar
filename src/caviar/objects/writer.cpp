
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
#include "caviar/CAVIAR.hpp"

#include "caviar/objects/writer.hpp"
#include "caviar/interpreter/communicator.hpp"

namespace caviar
{

    Writer::Writer(CAVIAR *fptr) : caviar_{fptr},
                                   comm{fptr->comm},
                                   error{fptr->error},
                                   output{fptr->output},
                                   input{fptr->input},
                                   object_handler{fptr->object_handler},
                                   object_container{fptr->object_container},
                                   object_creator{fptr->object_creator},
                                   log{fptr->log},
                                   in{fptr->in},
                                   out{fptr->out},
                                   err{fptr->err},
                                   log_flag{fptr->log_flag},
                                   out_flag{fptr->out_flag},
                                   err_flag{fptr->err_flag}

    {

        FC_OBJECT_INITIALIZE
        my_mpi_rank = fptr->comm->me;
        mpi_world_size = fptr->comm->nprocs;
    }

    Writer::~Writer()
    {
    }

    void Writer::verify_settings()
    {
    }

    void Writer::initialize() {}
    void Writer::write(int64_t, double) {}         // time_step and time
    void Writer::start_new_files() {}              // add_time_to_previous
    void Writer::start_new_files(std::string &) {} // add_time_to_previous
    void Writer::open_files() {}
    void Writer::close_files() {}
    void Writer::generate() {}

}
