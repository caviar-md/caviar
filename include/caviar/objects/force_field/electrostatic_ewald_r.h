
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATICEWALDR_H
#define CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATICEWALDR_H

#include "caviar/objects/force_field.h"

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  /**
   * This class electrostatic elwald real space force-field.
   *
   *
   */
  class Electrostatic_ewald_r : public Force_field
  {
  public:
    Electrostatic_ewald_r(class CAVIAR *);
    ~Electrostatic_ewald_r(){};

    double potential(const Vector<double> &);
    double potential(const int);

    Vector<double> field(const Vector<double> &);
    Vector<double> field(const int);

    bool read(class caviar::interpreter::Parser *);
    void verify_settings();
    void calculate_acceleration();

    double energy();

  public:
    double k_electrostatic, alpha;
    Vector<double> external_field;
  };

} // force_field

CAVIAR_NAMESPACE_CLOSE

#endif
