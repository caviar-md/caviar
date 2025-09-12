
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
#include "caviar/utility/vector.hpp"

namespace caviar
{

  namespace shape
  {
    namespace polyhedron
    {

      struct Polyhedron;
      /**
       * a utility class for polyhedrons
       */
      class Utility 
      {
      public:
        Utility(class CAVIAR *);
        virtual ~Utility();

        /**
         * after reading polyhedron file, it calculates normal vectors
         */
        void make_normal(shape::polyhedron::Polyhedron &);

        /**
         * makes normals of faces made of edges and other normals used in check_inside() algorithm.
         */
        void make_edge_norms(shape::polyhedron::Polyhedron &);

        /**
         * multiply all the normal Vectors with -1
         */
        void invert_normals(shape::polyhedron::Polyhedron &);

        /**
         * does what it says by using an inside point.
         */
        bool normals_are_pointing_outside(shape::polyhedron::Polyhedron &p_object, const Vector<double> &);
        FC_BASE_OBJECT_COMMON_TOOLS
      };
    } // polyhedron
  } // shape

}
