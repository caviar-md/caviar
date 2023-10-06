
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

#ifndef CAVIAR_OBJECTS_CONSTRAINT_NVE_H
#define CAVIAR_OBJECTS_CONSTRAINT_NVE_H

#include "caviar/objects/constraint.h"

CAVIAR_NAMESPACE_OPEN

namespace constraint
{

  /**
   * This class has N-V-E thermostat. It is done simply by velocity re-scaling.
   * It is the simplest way to control the energy in a MD simulation. However,
   * it won't gives a correct fluctuation for canonical ensembles.
   * The user has to set 'energy_tot' or 'energy_per_dof'.
   *
   * kinetic_energy (t)  = dof * k_b * T / 2   ,
   *  T : instantaneous temperature at time = t;
   *  dof : total degrees of freedom
   *
   * According to the formula above, we let the user to fix the temperatue if
   * wanted to do so. In that case, the users have two alternatives,
   * 1: set 'temperature' and 'k_b' (Boltzman constant)
   * 2: set 'k_b_t'
   *
   */
  class Nve : public Constraint
  {
  public:
    Nve(class CAVIAR *);
    ~Nve();
    bool read(class caviar::interpreter::Parser *);

    void apply_thermostat(int64_t, bool &recalculate_temperature);

    void verify_settings();

    double energy_per_dof, energy_tot;

    // if one has set_the temperature, it will use it a
    double temperature, kb, kbt;

    /**
     * Apply after each 'step' of timesteps passed
     */
    int step = 1;
  };

} // constraint

CAVIAR_NAMESPACE_CLOSE

#endif
