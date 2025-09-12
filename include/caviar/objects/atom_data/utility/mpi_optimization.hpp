
//========================================================================
//
// Copyright (C) 2023 by Morad Biagooi and Ehsan Nedaaee Oskoee.
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
#ifndef CAVIAR_ATOMDATA_MPIOPTIMIZATION_H
#define CAVIAR_ATOMDATA_MPIOPTIMIZATION_H
/**
 * More sharing results in less MPI overhead and probably faster MPI simulation. However, simulation will take more RAM. 
 * For large number of processors and particles, it must be tested.
*/
enum class MpiOptimization
{
  None,
  SingleMdDomain,
  ShareAtoms,
  ShareMolecules,
  ShareAtomsAndMolecules,
};

#endif