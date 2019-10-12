
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

#include "caviar/objects/unique/grid_1d.h"
#include "caviar/utility/interpreter_io_headers.h"

namespace caviar {
namespace objects {
namespace unique {


Grid_1D::Grid_1D (CAVIAR *fptr) : Unique{fptr},
    min{0}, max{0}, increment{-1},
    generated{false}, by_increment{false}, by_segment{false},
    segment{-1}, no_given_points{0} {
  FC_OBJECT_INITIALIZE_INFO
}
     
Grid_1D::Grid_1D (CAVIAR *fptr, double min, double max,
    double increment, int segment) : Unique{fptr},
    min{min}, max{max}, increment{increment},
    generated{false}, by_increment{false},
    by_segment{false}, segment{segment}, no_given_points{0} {
  FC_OBJECT_INITIALIZE_INFO
 }
  
Grid_1D::~Grid_1D () { }

void Grid_1D::verify_settings () {
  
}


bool Grid_1D::read (caviar::interpreter::Parser* parser) {
  FC_OBJECT_READ_INFO
    
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    if (token.string_value=="generate") {generate(); break;}      
    else ASSIGN_REAL(min,"GRID_1D Read: ","")
    else ASSIGN_REAL(max,"GRID_1D Read: ","")
    else ASSIGN_REAL(increment,"GRID_1D Read: ","")
    else ASSIGN_INT(segment,"GRID_1D Read: ","")  
    else error->all(FC_FILE_LINE_FUNC_PARSE,"Random_1D Read: Unknown variable or command ");
  }
  return in_file;;
}
  
void Grid_1D::generate () {
  output->info("Grid_1d Generate: ");    
  if (generated == true) 
    error->all("Grid_1D: Generate: cannot be generated twice. ");
  generated = true;
  if (segment>0 && increment>0)
    error->all("Grid_1D: Generate: Assigning both segment and increment is not possible. ");
  if (segment<0 && increment<0)
    error->all("Grid_1D: Generate: Assign one of segment or increment. ");
  if (min > max)
    error->all("Grid_1D: Generate: min has to be smaller than max. ");    
  if (segment<0) {
    by_increment = true;
    segment = int ((max-min)/increment);
  }
  if (increment<0) {
    by_segment = true;
    increment = (max - min)/double(segment);
  }
}
  
unsigned int Grid_1D::no_points () {
  if (by_segment)
    return segment + 1;
  else
    return segment;
}
  
double Grid_1D::give_point () {
  double val = min + no_given_points * increment;
  ++no_given_points;    
  if (by_segment) {
    if (no_given_points > segment) return max;
    else return val;
  } else {
    return val;    
  }
}
  
double Grid_1D::give_point (int i) {
  double val = min + i * increment;  
  if (by_segment) {
    //if (i == segment) return max; // XXX
    //else return val;
    return val;
  } else {
    return val;    
  }
}

} //unique
} //objects

} // namespace caviar

