
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

#include "caviar/objects/unique/random_1d.h"
#include "caviar/utility/interpreter_io_headers.h"

namespace caviar {

namespace unique {

//====================================
//====================================
//====================================

Random_1D::Random_1D (CAVIAR *fptr) : Unique{fptr},
   min{0}, max{2}, stddev{1}, mean{1}, type{"type"}, seed{1}, generated{false}
{
  FC_OBJECT_INITIALIZE_INFO
  type_int = 0;
}  
  
Random_1D::Random_1D (CAVIAR *fptr, std::string type, double min,
    double max, double stddev, double mean, int seed) : Unique{fptr},
    min{min}, max{max}, stddev{stddev}, mean{mean}, type{type}, seed{seed},
    generated{false}, generated_u_dist{false}, generated_n_dist{false}
{
  FC_OBJECT_INITIALIZE_INFO
  type_int = 0;
}
  
Random_1D::~Random_1D () {
  if (generated == true) 
    delete ran_gen;
  if (generated_u_dist == true)
    delete u_dist;
  if (generated_n_dist == true)
    delete n_dist;
}

void Random_1D::verify_settings () {
  
}


bool Random_1D::read (caviar::interpreter::Parser* parser) {
  FC_OBJECT_READ_INFO
    
    bool in_file = true;
    while(true) {
      GET_A_TOKEN_FOR_CREATION
      if (token.string_value=="GENERATE") {generate(); break;}
      else ASSIGN_STRING(type,"Random_1D creation: ","")
      else ASSIGN_REAL(min,"Random_1D creation: ","")
      else ASSIGN_REAL(max,"Random_1D creation: ","")
      else ASSIGN_REAL(stddev,"Random_1D creation: ","")
      else ASSIGN_REAL(mean,"Random_1D creation: ","")
      else ASSIGN_INT(seed,"Random_1D creation: ","")  
      else error->all(FC_FILE_LINE_FUNC_PARSE,"Random_1D creation: Unknown variable or command ");
    }
    return in_file;;

  }
  
void Random_1D::generate () {
  output->info("Random_1D Generate: ");  

  if (generated == true) 
    error->all(FC_FILE_LINE_FUNC, "Generate: cannot be generated twice. ");
  generated = true;  

  ran_gen = new std::mt19937 (seed);

  if (type=="UNIFORM") {
    type_int = 1;
    u_dist = new std::uniform_real_distribution<> (min, max) ;
    generated_u_dist = true;
  } else if (type=="NORMAL") {
    type_int = 2;
    n_dist = new std::normal_distribution<> (mean, stddev);
    generated_n_dist = true;
  } else error->all (FC_FILE_LINE_FUNC, "Expected NORMAL or UNIFORM for the type. ");  
}

double Random_1D::give_value () {
  if (!generated)
    error->all(FC_FILE_LINE_FUNC, "This class is not generated yet. ");
  error->all(FC_FILE_LINE_FUNC, "not implemented. ");
  switch(type_int) {
    case (1):
//      return *u_dist(*ran_gen);

    case (2):
//      return *n_dist(*ran_gen);

    default:
      error->all(FC_FILE_LINE_FUNC,"undefined random distribution.");
  }
  return 1.0;
}
  
} //unique


} // namespace caviar

