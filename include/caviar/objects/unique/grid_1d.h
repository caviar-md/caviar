
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

#ifndef CAVIAR_OBJECTS_UNIQUE_GRID1D_H
#define CAVIAR_OBJECTS_UNIQUE_GRID1D_H

#include "caviar/objects/unique.h"

CAVIAR_NAMESPACE_OPEN
class Parser;

namespace unique {

/**
 * This class creates grid positions for the initial position of the particles.
 * 
 */
class Grid_1D : public Unique {
  public:
    Grid_1D (class CAVIAR *) ;    
    Grid_1D (class CAVIAR *, double MIN, double MAX, double increment, int segment) ;
    ~Grid_1D () ;
    void verify_settings ();    
    bool read (caviar::interpreter::Parser *);
    void generate (); // calculates the parameters
    unsigned int no_points ();
    double give_point ();
    double give_point (int);
    
    void reset();

    double min, max, increment;

    bool generated; // true if generate() has been called.    
    bool by_increment, by_segment;
    
    int segment;    
    int no_given_points;     
    int num; // number of random atoms or molecules to be created        
    int type_int;
   
    std::string TYPE;
    
};

} //unique


CAVIAR_NAMESPACE_CLOSE

#endif
