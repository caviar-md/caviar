
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Implementation of Shake algorithm is done by Ashkan Shahmoradi
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

#ifndef CAVIAR_OBJECTS_CONSTRAINT_SHAKE_H
#define CAVIAR_OBJECTS_CONSTRAINT_SHAKE_H

#include "caviar/objects/constraint.h"

CAVIAR_NAMESPACE_OPEN

class Domain;
namespace constraint
{

  /**
   * This class fixes atomic bonds using SHAKE method
   *
   *
   */
  class Shake : public Constraint
  {
  public:
    Shake(class CAVIAR *);
    ~Shake();
    bool read(class caviar::interpreter::Parser *);

    void apply_shake(int64_t);

    void update_velocity_after_shake(int64_t);

    void bond_fix();

    void verify_settings();

    static inline int delta(int a, int b)
    {
      if (a == b)
        return 1;
      else
        return 0;
    }

    class Domain *domain = nullptr;

    /**
     * The molecule bond types which will be included in the algorithm
    */    
    std::vector<int> bond_type;

    double dt;
    double error_tolerance;

    bool initialized;

    int iteration_max = 100;
  };

} // constraint

CAVIAR_NAMESPACE_CLOSE

#endif
