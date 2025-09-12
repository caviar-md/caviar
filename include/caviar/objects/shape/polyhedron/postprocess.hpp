
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

#pragma once

#include <string>
#include "caviar/utility/objects_common_headers.hpp"
#include "caviar/objects/shape/polyhedron.hpp"

namespace caviar
{

  namespace shape
  {
    namespace polyhedron
    {
      struct Polyhedron;
      class Postprocess 
      {
      public:
        Postprocess(class CAVIAR *);
        virtual ~Postprocess();

        /**
         * contains the faces neccesary to check
         * make_grid has to be called after lowest_highest_coord()
         */
        void make_grid(shape::polyhedron::Polyhedron &);

        /**
         * calculates gxlo, gxhi, gylo...
         */
        void lowest_highest_coord(shape::polyhedron::Polyhedron &);
        FC_BASE_OBJECT_COMMON_TOOLS
      };
    } // polyhedron
  } // shape

}
