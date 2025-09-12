
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

#include "caviar/objects/shape/polyhedron/input.hpp"
#include "caviar/objects/shape/polyhedron/format_unv_reader.hpp"
#include "caviar/objects/shape/polyhedron/format_vtk_reader.hpp"
#include "caviar/objects/shape/polyhedron/format_stl_reader.hpp"
#include "caviar/CAVIAR.hpp"

#include <string>
#include <fstream>

namespace caviar
{

  namespace shape
  {
    namespace polyhedron
    {

      Input::Input(CAVIAR *fptr) : caviar_{fptr},
                                   comm{fptr->comm},
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

      Input::~Input() {}
      void Input::verify_settings() {}
      void Input::read_vtk(shape::polyhedron::Polyhedron &p_object, const std::string &file_name)
      {
        class Format_vtk_reader fvr(caviar_);
        fvr.read_polyhedron(p_object, file_name);
      }

      void Input::read_stl(shape::polyhedron::Polyhedron &p_object, const std::string &file_name)
      {
        class Format_stl_reader fvr(caviar_);
        fvr.read_polyhedron(p_object, file_name);
      }

      void Input::read_unv(shape::polyhedron::Polyhedron &p_object, const std::string &file_name)
      {
        class Format_unv_reader fvr(caviar_);
        fvr.read_polyhedron(p_object, file_name);
      }

    } // polyhedron
  } // shape

}
