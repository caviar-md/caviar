
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

#include "point_condition.h"
#include "macros.h"
#include <iostream>
//#include <cmath>

#define DEFAULT_TOL 1E-5
#define ABS(x) (x<0) ? -x : x 

Point_condition::Point_condition() {
  std::cout << ERROR_COMMENT << "is not supported.\n";
}

Point_condition::~Point_condition() {

}

Point_condition::Point_condition(const std::string el, const std::string op, double v) {
  if (el == "x") {
    vector_component = 0;
  } else   if (el == "y") {
    vector_component = 1;
  } else   if (el == "z") {
    vector_component = 2;
  } else  {
      std::cout << ERROR_COMMENT << "'vector_component'" << el << " is not supported\n";      
    return;    
  } 

  set_operator(op, operation_type);

  value = v;
  condition_type = 0;
}

Point_condition::Point_condition(const std::string st, const Vector<double> vec_, const std::string op, double v) {
  if (st == "distance") {
    condition_type = 1;
  } else  {
      std::cout << ERROR_COMMENT << "command " << st << " is not supported\n";      
    return;    
  } 

  set_operator(op, operation_type);

  value = v;
  value_sq = value*value;
  vec = vec_;
}

void Point_condition::set_operator(const std::string op, int &type) {
  if (op=="<" || op=="smaller") {
    type = -2;
  } else   if (op=="<=" || op=="eqsmaller") {
    type = -1;
  } else  if (op=="==" || op=="equal") {
    type = 0;
  } else  if (op==">=" || op=="eqlarger") {
    type = 1;
  } else  if (op==">" || op=="larger") {
    type = 2;
  } else  {
        std::cout << ERROR_COMMENT << "'operation_type'" << op << " is not supported\n";
    return;
  }
}

bool Point_condition::in_condition(const Vector<double> &v) const {
  
  switch (condition_type) {

  // comparing one condition
  case 0: {
   
    switch (vector_component) {
    case 0:{
      switch (operation_type) {
      case -2:
        return (v.x <  value);

      case -1:
        return (v.x <= value);

      case 0:
        //return (v.x == value);
        return (ABS(v.x - value) < DEFAULT_TOL );

      case +1:
        return (v.x >= value);

      case +2:
        return (v.x >  value);

      default:
        std::cout << ERROR_COMMENT << "'operation_type'" << operation_type << " is not supported\n";
        return false;
      }
    } 

    case 1:{
      switch (operation_type) {
      case -2:
        return (v.y <  value);

      case -1:
        return (v.y <= value);

      case 0:
        //return (v.y == value);
        return (ABS(v.y - value) < DEFAULT_TOL );

      case +1:
        return (v.y >= value);

      case +2:
        return (v.y >  value);

      default:
        std::cout << ERROR_COMMENT << "'operation_type'" << operation_type << " is not supported\n";
        return false;
      }      
    } 


    case 2:{
      switch (operation_type) {
      case -2:
        return (v.z <  value);

      case -1:
        return (v.z <= value);

      case 0:
        //return (v.z == value);
        return (ABS(v.z - value) < DEFAULT_TOL );

      case +1:
        return (v.z >= value);

      case +2:
        return (v.z >  value);

      default:
        std::cout << ERROR_COMMENT << "'operation_type'" << operation_type << " is not supported\n";
        return false;
      }      
    } 

    default:{
      std::cout << ERROR_COMMENT << "'vector_component'" << vector_component << " is not supported\n";      
      return false;
    } 
    }

  } 


  case 1:{
    auto v_dif = v - vec;
    auto v_dif_sq = v_dif*v_dif;
    switch (operation_type) {
    case -2:
      return (v_dif_sq <  value_sq);

    case -1:
      return (v_dif_sq <= value_sq);

    case 0:
      //return (v_dif_sq == value_sq);
      return (ABS(v_dif_sq - value_sq) < DEFAULT_TOL );

    case +1:
      return (v_dif_sq >= value_sq);

    case +2:
      return (v_dif_sq >  value_sq);

    default:
      std::cout << ERROR_COMMENT << "'operation_type'" << operation_type << " is not supported\n";
      return false;
    }
  } 


  default:{
    std::cout << ERROR_COMMENT << "'condition_type'" << condition_type << " is not supported\n";
    return false;
  } 

  }

}
