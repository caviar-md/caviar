
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

#include "caviar/objects/unique/time_function.h"
#include "caviar/utility/interpreter_io_headers.h"


#if defined(CAVIAR_WITH_MUPARSER)
  #include <muParser.h>
#endif
CAVIAR_NAMESPACE_OPEN

namespace unique {

//====================================
//====================================
//====================================

Time_function::Time_function (CAVIAR *fptr) : Unique{fptr}
{
  FC_OBJECT_INITIALIZE_INFO
  export_values_to_file = false;
  export_file_append = false;
#if defined(CAVIAR_WITH_MUPARSER)  
  muParser = new mu::Parser;
  muParser->DefineVar("t", &time_variable); 
#else
  error->all(FC_FILE_LINE_FUNC,"CAVIAR must be build with 'CAVIAR_WITH_MUPARSER=ON' for 'Time_function'");
#endif
}  
  
  
Time_function::~Time_function () {
#if defined(CAVIAR_WITH_MUPARSER)
  delete muParser;

#endif  
  if (export_values_to_file)
    if (ofs_time_value.is_open())    ofs_time_value.close();  
}

void Time_function::verify_settings () {
  
}


bool Time_function::read (caviar::interpreter::Parser* parser) 
{
  FC_OBJECT_READ_INFO
    
    bool in_file = true;
    while(true) {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;      
      if (string_cmp(t,"set_formula")) {
        const auto t2 = parser->get_val_token();
        function_definition = t2.string_value;
        generate_formula();
        continue;
      } else if (string_cmp(t,"export_file_name")) {
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
      else error->all(FC_FILE_LINE_FUNC_PARSE,"Time_function creation: Unknown variable or command ");
    }
    return in_file;;

}
 
  
void Time_function::generate_formula () {	  
#if defined(CAVIAR_WITH_MUPARSER)
  muParser->SetExpr(function_definition);                                         
#endif
}

void Time_function::generate_export_file () {	  
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

void Time_function::update_time_variable(double t) {

  time_variable = t;
  calculate ();
  if (export_values_to_file)
  {
    ofs_time_value << t << " " << current_value << "\n";
  }  
}
  
void Time_function::calculate ()
{
#if defined(CAVIAR_WITH_MUPARSER)
  current_value = muParser->Eval();  
#endif
}
  
} //unique


CAVIAR_NAMESPACE_CLOSE

