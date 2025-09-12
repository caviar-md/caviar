
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
#include "caviar/CAVIAR.hpp"
#include "caviar/interpreter/all.hpp"
#include "caviar/interpreter/object_container.hpp"
#include "caviar/utility/interpreter_io_headers.hpp"
#include "caviar/objects/all.hpp" // used in deletion

namespace caviar
{
  namespace interpreter
  {

    Object_container::Object_container(CAVIAR *fptr) : caviar_{fptr}, comm{fptr->comm},
                                                       error{fptr->error},
                                                       output{fptr->output},
                                                       input{fptr->input},
                                                       object_handler{fptr->object_handler},
                                                       object_container{fptr->object_container},
                                                       object_creator{fptr->object_creator},
                                                       log{fptr->log},
                                                       in{fptr->in},
                                                       out{fptr->out},
                                                       err{fptr->err},
                                                       log_flag{fptr->log_flag},
                                                       out_flag{fptr->out_flag},
                                                       err_flag{fptr->err_flag} {}

    Object_container::~Object_container()
    {

#define FC_GENERAL_CLASSNAME_MACRO(VAR1, VAR2, VAR3) \
  for (unsigned int i = 0; i < VAR2.size(); i++)     \
  {                                                  \
    delete (VAR2[i]);                                \
  }                                                  \
  VAR2.clear();

#define FC_GENERAL_CLASSNAME_MACRO_ACTIVATED

#include "caviar/objects/macro/declaration/all.hpp"

#undef FC_GENERAL_CLASSNAME_MACRO_ACTIVATED
#undef FC_GENERAL_CLASSNAME_MACRO
    }
    void Object_container::verify_settings() {}
    bool Object_container::read(caviar::interpreter::Parser *parser)
    {
      output->info("object_container read");
      bool in_file = true;

      while (true)
      {
        GET_A_TOKEN_FOR_CREATION
        auto t = token.string_value;
        if (string_cmp(t, "report"))
        {
          report();
        }
        else
          error->all(FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
      }

      return in_file;
    }

    void Object_container::report()
    {
      output->info("object_container report:");
      std::unordered_set<std::string>::iterator it;
      for (it = all_names.begin(); it != all_names.end(); ++it)
        std::cout << *it << std::endl;
    }

  } // interpreter
}
