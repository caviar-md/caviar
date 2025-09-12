
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

#ifndef CAVIAR_OBJECTS_UNIQUE_RANDOM1D_H
#define CAVIAR_OBJECTS_UNIQUE_RANDOM1D_H

#include "caviar/objects/unique.h"
#include <random>

CAVIAR_NAMESPACE_OPEN
class Parser;

namespace unique
{

  /**
   * This class is a wrapper for std::random used for initial position of the particles
   *
   */
  class Random_1D : public Unique
  {
  public:
    Random_1D();
    Random_1D(class CAVIAR *);
    Random_1D(class CAVIAR *, std::string TYPE, double MIN, double MAX, double STDDEV, double MEAN, int SEED);
    ~Random_1D();
    bool read(caviar::interpreter::Parser *);
    void generate(); // creates the std::mt19937 with given parameters ...
    void verify_settings();
    double give_value();
    double min, max, stddev, mean;

    std::string type;

    int seed;
    int num; // number of random atoms or molecules to be created
    int type_int;

    bool generated, generated_u_dist, generated_n_dist;

    std::mt19937 *ran_gen;
    std::uniform_real_distribution<> *u_dist;
    std::normal_distribution<> *n_dist;
  };

} // unique

CAVIAR_NAMESPACE_CLOSE

#endif
