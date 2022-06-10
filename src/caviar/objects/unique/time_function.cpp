
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

#include <muParser.h>

namespace caviar {
namespace objects {
namespace unique {

//====================================
//====================================
//====================================

Time_function::Time_function (CAVIAR *fptr) : Unique{fptr}
{
  FC_OBJECT_INITIALIZE_INFO
  muParser = new mu::Parser;
  muParser->DefineVar("t", &time_variable);         
}  
  
Time_function::Time_function (CAVIAR *fptr, std::string f) : Unique{fptr},
    function_definition{f}
{
  FC_OBJECT_INITIALIZE_INFO
  muParser = new mu::Parser;
  muParser->DefineVar("t", &time_variable);           
}
  
Time_function::~Time_function () {
  delete muParser;
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
      } else if (string_cmp(t,"export_file_name")) {
        const auto t2 = parser->get_val_token();
        export_file_name = token.string_value;
      } 
      else error->all(FC_FILE_LINE_FUNC_PARSE,"Time_function creation: Unknown variable or command ");
    }
    return in_file;;

}
  
void Time_function::generate () {	  
  muParser->SetExpr(function_definition);                                         
}

void Time_function::calculate ()
{
  muParser->Eval();  
}
  
} //unique
} //objects

} // namespace caviar

