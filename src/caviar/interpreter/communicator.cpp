
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

#include "caviar/interpreter/communicator.h"
#include "caviar/utility/vector.h"
#if defined(CAVIAR_WITH_MPI)
#include <mpi.h>
#endif

CAVIAR_NAMESPACE_OPEN
namespace interpreter
{
  Communicator::Communicator(CAVIAR *fptr) : Pointers{fptr}
  {
#if defined(CAVIAR_WITH_MPI)

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /*
    // Adding a MPI type for caviar::Vector
    const int nitems = 3;
    int blocklengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    // MPI_Datatype mpi_fc_vector_type;
    MPI_Aint offsets[3];

    offsets[0] = offsetof(caviar::Vector<double>, x);
    offsets[1] = offsetof(caviar::Vector<double>, y);
    offsets[2] = offsetof(caviar::Vector<double>, z);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_fc_vector_type);
    MPI_Type_commit(&mpi_fc_vector_type);
    */
#endif
  }

  void Communicator::broadcast(bool &flag)
  {
#if defined(CAVIAR_WITH_MPI)
    MPI_Bcast(&flag, 1, MPI::BOOL, 0, MPI_COMM_WORLD);
#else
    std::cout << "Communicator::broadcast " << flag << std::endl;
#endif
  }

  void Communicator::broadcast(size_t &n)
  {
#if defined(CAVIAR_WITH_MPI)
    MPI_Bcast(&n, 1, MPI::INT, 0, MPI_COMM_WORLD);
#else
    std::cout << "Communicator::broadcast " << n << std::endl;
#endif
  }

  void Communicator::broadcast(size_t &n, char *str)
  {
#if defined(CAVIAR_WITH_MPI)
    MPI_Bcast(&n, 1, MPI::INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(str, n, MPI::CHAR, 0, MPI_COMM_WORLD);
#else
    std::cout << "Communicator::broadcast " << n << " " << str << std::endl;
#endif
  }

  void Communicator::broadcast(std::string &str)
  {
#if defined(CAVIAR_WITH_MPI)

    int n = me == 0 ? str.length() : 0;
    MPI_Bcast(&n, 1, MPI::INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // The reason behind '[n+1]' for 'char *tmp' is because of 'strcpy':
    // (from http://www.cplusplus.com/reference/cstring/strcpy/)
    // To avoid overflows, the size of the array pointed by destination shall be
    // long enough to contain the same C string as source (including the
    // terminating null character), and should not overlap in memory with source.
    char *tmp = new char[n + 1];

    strcpy(tmp, str.c_str());

    MPI_Bcast(tmp, n, MPI::CHAR, 0, MPI_COMM_WORLD);

    if (tmp)
    {
      if (tmp[0])
      {
        str.assign(const_cast<const char *>(tmp), n);
      }
      else
      {
      }
    }
    else
    {
    }

    MPI_Barrier(MPI_COMM_WORLD);

    delete[] tmp;

#else
    std::cout << "Communicator::broadcast is called in non-mpi mode" << str << std::endl;
#endif
  }

} // interpreter
CAVIAR_NAMESPACE_CLOSE
