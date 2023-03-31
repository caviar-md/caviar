
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

#include "caviar/objects/shape/mixed.h"
#include "caviar/utility/interpreter_io_headers.h"

namespace caviar {

namespace shape {

Mixed::Mixed (CAVIAR *fptr) : Shape {fptr} {
  FC_OBJECT_INITIALIZE_INFO
}
   
Mixed::~Mixed () {}  
  
bool Mixed::read (caviar::interpreter::Parser* parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

#define AND_OR_INSIDE_OUTSIDE \
if (shapes.size() == 0) {\
  if (string_cmp(t,"inside")) inside_check = true;\
  else if (string_cmp(t,"outside")) inside_check = false;\
  else error->all(FC_FILE_LINE_FUNC_PARSE,"expected 'INSIDE' or 'OUTSIDE': ");\
} else if (string_cmp(t,"inside") || string_cmp(t,"outside"))\
    error->all(FC_FILE_LINE_FUNC_PARSE,"expected 'AND_INSIDE','AND_OUTSIDE', 'OR_INSIDE','OR_OUTSIDE': ");\
std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator it_1;\
std::string name_1;\
GET_A_STRING(name_1,""," expected an Mixed or SHAPE NAME: ")\
CHECK_NAME_EXISTANCE(name_1, it_1, "","")\
if (it_1->second.type == caviar::interpreter::object_handler::gdst("shape"))\
  shapes.push_back(object_container->shape[it_1->second.index]);\
else error->all(FC_FILE_LINE_FUNC_PARSE,"expected a SHAPE NAME.");
          
  while(true) {
    
      GET_A_TOKEN_FOR_CREATION
      const auto t = token.string_value;
      if (string_cmp(t,"inside")) {        
        AND_OR_INSIDE_OUTSIDE
        operators.push_back (0);  //inside_check = true     
      } else if (string_cmp(t,"outside")) {        
        AND_OR_INSIDE_OUTSIDE
        operators.push_back (0); //inside_check = flase           
      } else if (string_cmp(t,"and_inside")) {        
        AND_OR_INSIDE_OUTSIDE
        operators.push_back (1);       
      } else if (string_cmp(t,"and_outside")) { 
        AND_OR_INSIDE_OUTSIDE
        operators.push_back (-1);
      }  else if (string_cmp(t,"or_inside")) { 
        AND_OR_INSIDE_OUTSIDE
        operators.push_back (2);
      }  else if (string_cmp(t,"or_outside")) { 
        AND_OR_INSIDE_OUTSIDE
        operators.push_back (-2);
      }              
      else error->all(FC_FILE_LINE_FUNC_PARSE,"Mixed Read: Unknown variable or command ");
    }
    return in_file;;
#undef AND_OR_INSIDE_OUTSIDE
}
  
bool Mixed::is_inside (const Vector<double> &v) {
    bool tmp;
    
    if (inside_check) tmp = shapes[0]->is_inside(v);
    else tmp = !shapes[0]->is_inside(v);
    
    for (unsigned int i = 1; i < operators.size(); ++i) {
      switch (operators[i]) {
        case 1:
          tmp = (shapes[i]->is_inside(v) && tmp);
        break;

        case -1:
          tmp = (shapes[i]->is_outside(v) && tmp);                
        break;
        
        case 2:
          tmp = (shapes[i]->is_inside(v) || tmp);                
        break;
        
        case -2:
          tmp = (shapes[i]->is_outside(v) || tmp);                
        break;
        
        default:
          //error->all(FC_FILE_LINE_FUNC_PARSE,"Mixed is_inside: undefined boolean operator. ");        XXX
        break;                
      }
    }
    return tmp;
}
  
bool Mixed::is_inside (const Vector<double> &v, const double r) {
    bool tmp;
    
    if (inside_check) tmp = shapes[0]->is_inside(v, r);
    else tmp = !shapes[0]->is_inside(v, r);
    
    for (unsigned int i = 1; i < operators.size(); ++i) {
      switch (operators[i]) {
        case 1:
          tmp = (shapes[i]->is_inside(v, r) && tmp);
        break;

        case -1:
          tmp = (shapes[i]->is_outside(v, r) && tmp);                
        break;
        
        case 2:
          tmp = (shapes[i]->is_inside(v, r) || tmp);                
        break;
        
        case -2:
          tmp = (shapes[i]->is_outside(v, r) || tmp);                
        break;
        
        default:
          //error->all(FC_FILE_LINE_FUNC_PARSE,"Mixed is_inside: undefined boolean operator. ");    XXX    
        break;                
      }
    }
    return tmp;
}
  
  


bool Mixed::in_contact(const Vector<double> &v, const double r, Vector<double> & contact_vector) {
  std::string s = "incomplete function:";
  s += __FILE__ + std::to_string(__LINE__) + __func__;  
  output->warning(s);
  std::cout << "  " << v << r << contact_vector   << std::endl;     
  return false;
}

} //shape

} // namespace caviar

