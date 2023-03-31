
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

#include "caviar/objects/unique/time_function_3d.h"
#include "caviar/utility/interpreter_io_headers.h"


#if defined(CAVIAR_WITH_MUPARSER)
  #include <muParser.h>
#endif
CAVIAR_NAMESPACE_OPEN

namespace unique {

//====================================
//====================================
//====================================

Time_function_3d::Time_function_3d (CAVIAR *fptr) : Unique{fptr}
{
  FC_OBJECT_INITIALIZE_INFO
  export_values_to_file = false;
  export_file_append = false;
#if defined(CAVIAR_WITH_MUPARSER)  
  muParser_x = new mu::Parser;
  muParser_y = new mu::Parser;
  muParser_z = new mu::Parser;
  muParser_x->DefineVar("t", &time_variable); 
  muParser_y->DefineVar("t", &time_variable); 
  muParser_z->DefineVar("t", &time_variable); 

#else
  error->all(FC_FILE_LINE_FUNC,"CAVIAR must be build with 'CAVIAR_WITH_MUPARSER=ON' for 'Time_function_3d'");
#endif
}  
  
  
Time_function_3d::~Time_function_3d () {
#if defined(CAVIAR_WITH_MUPARSER)
  delete muParser_x;
  delete muParser_y;
  delete muParser_z;
#endif  
  if (export_values_to_file)
    if (ofs_time_value.is_open())    ofs_time_value.close();  
}

void Time_function_3d::verify_settings () {
  
}


bool Time_function_3d::read (caviar::interpreter::Parser* parser) 
{
  FC_OBJECT_READ_INFO
    
    bool in_file = true;
    while(true) {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;      
      if (string_cmp(t,"set_formula_x")) {
        const auto t2 = parser->get_val_token();
        function_definition_x = t2.string_value;
        generate_formula();
        continue;
      } 
      else if (string_cmp(t,"set_formula_y")) {
        const auto t2 = parser->get_val_token();
        function_definition_y = t2.string_value;
        generate_formula();
        continue;
      }
      else if (string_cmp(t,"set_formula_z")) {
        const auto t2 = parser->get_val_token();
        function_definition_z = t2.string_value;
        generate_formula();
        continue;
      }else if (string_cmp(t,"export_file_name")) {
        const auto t2 = parser->get_val_token();
        export_file_name = t2.string_value;
        if (export_file_name=="")
          export_values_to_file = false;        
        else
          export_values_to_file = true;
        generate_export_file ();
        continue;
      } else if (string_cmp(t,"export_file_append")) {
        export_file_append = true;
      }
      else error->all(FC_FILE_LINE_FUNC_PARSE,"Time_function_3d creation: Unknown variable or command ");
    }
    return in_file;;

}
 
  
void Time_function_3d::generate_formula () {	  
#if defined(CAVIAR_WITH_MUPARSER)
  muParser_x->SetExpr(function_definition_x);  
  muParser_y->SetExpr(function_definition_y);                                         
  muParser_z->SetExpr(function_definition_z);                                                                                
#endif
}

void Time_function_3d::generate_export_file () {	  
  // if script changes the parameter to "" in middle of simulation
  // the file will be closed
  if (!export_values_to_file)
  {
    if (ofs_time_value.is_open())
      ofs_time_value.close();
  }
  
  // opening/appending a file to export
  if (export_values_to_file)
  {     
    if (!ofs_time_value.is_open()) 
    {
      if (export_file_append)
        ofs_time_value.open(export_file_name.c_str(), std::ios_base::app);
      else
        ofs_time_value.open(export_file_name.c_str());
    }
  }
}

void Time_function_3d::update_time_variable(double t) {

  time_variable = t;
  calculate ();
  if (export_values_to_file)
  {
    ofs_time_value << t << " " << current_value.x << " " << current_value.y << " " << current_value.z << "\n";
  }  
}
  
void Time_function_3d::calculate ()
{
#if defined(CAVIAR_WITH_MUPARSER)
  current_value.x = muParser_x->Eval();  
  current_value.y = muParser_y->Eval();  
  current_value.z = muParser_z->Eval(); 
#endif

}
  
} //unique


CAVIAR_NAMESPACE_CLOSE

