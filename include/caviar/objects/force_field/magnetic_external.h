
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_MAGNETICEXTERNAL_H
#define CAVIAR_OBJECTS_FORCEFIELD_MAGNETICEXTERNAL_H

#include "caviar/objects/force_field.h"

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  /**
   * This class calculates a very simple external magnetic force-field for the charged
   *  particles with velocity as 'F = q V x B'
   */
  class Magnetic_external : public Force_field
  {
  public:
    Magnetic_external(class CAVIAR *);
    ~Magnetic_external(){};

    bool read(class caviar::interpreter::Parser *);
    void verify_settings();
    void calculate_acceleration();

  public:
    double amplitude;
    Vector<double> direction;
  };

} // force_field

CAVIAR_NAMESPACE_CLOSE

#endif
