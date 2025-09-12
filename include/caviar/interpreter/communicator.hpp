
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

#pragma once

#include "caviar/utility/interpreter_common_headers.hpp"
// #if defined(CAVIAR_WITH_MPI)
// #include<mpi.h>
// #endif

#if defined(CAVIAR_WITH_MPI)
   // MPI_Datatype mpi_fc_vector_type;
#endif

namespace caviar
{
namespace interpreter
{

  /**
   * This class handles MPI process base communications.
   *
   */
  class Communicator 
  {
  public:
    Communicator(CAVIAR *);
    virtual ~Communicator();

    // broadcast a variable from the root process to others
    void broadcast(bool &);
    void broadcast(size_t &);
    void broadcast(size_t &, char *);
    void broadcast(std::string &);

#if defined(CAVIAR_WITH_MPI)
    //MPI_Datatype mpi_fc_vector_type;
#endif

    int me, nprocs; // MPI process rank and number of processes
    FC_BASE_OBJECT_COMMON_TOOLS
  };
} // interpreter
}
