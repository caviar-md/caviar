
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
#include <fstream>

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

        void apply_barostat(int64_t, bool &fix_position_needed);

        void verify_settings();

        double get_pressure();
        
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
        int step = 5;

        /**
         * Moving Average window
        */
        int ma_window = 5;

        /**
         * Moving Average Type. 0: No moving average, 1: Simple moving average
        */
        int ma_type = 0;


        /**
         * export data for parameters validation
        */
        bool export_data = false;


        /**
         * export data for parameters validation
        */
        std::ofstream ofs_export;
        
        /**
         * Moving Average counter
        */
        int ma_counter = 1;

        /**
         * Moving Average current value
        */
        double ma_current = 0.0;


        

        /**
         * Maximum value of scaling
        */
        double xi_max = 0.1;


        /**
         * The direction of scaling simulation box. It will be done on periodic direction of the domain.
        */
        caviar::Vector<int> scale_axis;


        /**
        * Function to update moving average
        */
        double updateMovingAverage(double newData, double prevAvg, int numDataPoints); 


    public:
    };

} // constraint

CAVIAR_NAMESPACE_CLOSE

#endif
