
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

#include "caviar/utility/objects_common_headers.hpp"
#include "caviar/objects/shape/polyhedron.hpp"

namespace caviar
{

  namespace shape
  {
    namespace polyhedron
    {
      struct Polyhedron;
      class Format_vtk_reader 
      {
      public:
        Format_vtk_reader(class CAVIAR *);
        virtual ~Format_vtk_reader();

        void read_polyhedron(shape::polyhedron::Polyhedron &, const std::string &);
        void write_unstructured_vtk4(shape::polyhedron::Polyhedron &, const std::string st_out = "o_shape_unstructured.vtk");
        void write_polydata_vtk4(shape::polyhedron::Polyhedron &, const std::string st_out = "o_shape_polydata.vtk");

        // It checks if the two vertices are similar
        // Then makes a map of all vertices to the similar
        // ones with the lower index, then clear
        // void merge_vertices (int);
        FC_BASE_OBJECT_COMMON_TOOLS
      };
    } // polyhedron
  } // shape

}
