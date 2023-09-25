
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

#ifndef CAVIAR_OBJECTS_CONSTRAINT_BERENDSENBAROSTAT_H
#define CAVIAR_OBJECTS_CONSTRAINT_BERENDSENBAROSTAT_H

#include "caviar/objects/constraint.h"

CAVIAR_NAMESPACE_OPEN

class Force_field;
class Domain;

namespace constraint
{
    
    /**
     * Brendsen thermostat implementation
     * Xi = 1 - kappa (dt / tp) (P - P_int(t)) 
     * r' = Xi^(1/3) * r
     * for liquid water: tp = 0.01 ps  ...  tp = 0.1 ps
     * for liquid water: kappa ~ 4.5 * 10E-5 (Bar^(-1))
     */
    class Berendsen_barostat : public Constraint
    {
    public:
        Berendsen_barostat(class CAVIAR *);
        ~Berendsen_barostat();
        bool read(class caviar::interpreter::Parser *);

        void apply_on_velocity(int64_t);

        void verify_settings();
        
        /**
         * the force_fields which a scale_position() must be called.
         * Not all force_fields needs that. 
        */        
        std::vector<Force_field *> force_field;

        Domain *domain = nullptr;

        double tp;
        double dt;
        double kappa;
        double pressure;

        /**
         * Apply after each 'step' of timesteps passed
        */
        int step = 1;

        /**
         * Maximum value of scaling
        */
        double xi_max = 0.0001;


        /**
         * The direction of scaling simulation box. It will be done on periodic direction of the domain.
        */
        caviar::Vector<int> scale_axis;

    public:
    };

} // constraint

CAVIAR_NAMESPACE_CLOSE

#endif
