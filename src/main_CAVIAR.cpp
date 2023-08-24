
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

#ifdef CAVIAR_WITH_MPI
#include <mpi.h>
#endif

#include "caviar/CAVIAR.h"

#ifdef CAVIAR_WITH_DEALII_MPI
#include <deal.II/base/mpi.h>
#endif

int main(int argc, char **argv)
{

#if defined(CAVIAR_WITH_DEALII_MPI)
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
#elif defined(CAVIAR_WITH_MPI)
  MPI_Init(&argc, &argv);
#endif

  caviar::CAVIAR caviar(argc, argv);

  try
  {
    caviar.execute();
  }
  catch (std::exception &exc)
  {
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << "Unknown exception!" << std::endl;
  }

#if defined(CAVIAR_WITH_DEALII_MPI)
  // Do nothing
#elif defined(CAVIAR_WITH_MPI)
  MPI_Finalize();
#endif

  return 0;
}
