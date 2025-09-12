
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


#include "caviar/utility/vector3d.hpp"
#include "caviar/objects/shape/polyhedron/polyhedron.hpp"
#include "caviar/utility/objects_common_headers.hpp"

#include <vector>

namespace caviar
{

  namespace shape
  {
    namespace polyhedron
    {
      class Input;
      class Preprocess;
      class Postprocess;
      class Utility;
      class Point_Inside;
      class Output;
      class Handler 
      {
      public:
        Handler(class CAVIAR *);
        virtual ~Handler();

        bool read(caviar::interpreter::Parser *);
        bool is_inside(const Vector3d<double> &v);
        bool is_inside(const Vector3d<double> &, const double rad);
        bool in_contact(const Vector3d<double> &, const double rad, Vector3d<double> &contact_vector);

      public:
        void command_generate();

        class shape::polyhedron::Input *polyhedron_input;
        class shape::polyhedron::Preprocess *polyhedron_preprocess;
        class shape::polyhedron::Postprocess *polyhedron_postprocess;
        class shape::polyhedron::Utility *polyhedron_utility;
        class shape::polyhedron::Point_Inside *polyhedron_point_inside;
        class shape::polyhedron::Output *polyhedron_output;

        struct caviar::shape::polyhedron::Polyhedron polyhedron;

        // this variable will be set for polyhedron.grid_tol .
        double radius_max;

        bool polyhedron_read, output_mesh_tcl, output_normals_tcl, output_edges_tcl,
            output_mesh_povray, output_normals_vectors;
        bool invert_normals, correct_normals, use_grid;

        /**
         * the class uses this point to decide if the overall direction of the normals
         * are pointing inside. The current algorithm has a weakness.
         * that is if the geometry is complex, (if it is concave
         * near a corner) the user has to set the point to a place which
         * is near to the center of a flat polyhedron, and far from concave corners.
         */
        Vector3d<double> an_inside_point;
        bool an_inside_point_is_set;
        FC_BASE_OBJECT_COMMON_TOOLS
      };
    } // polyhedron
  } // shape

}
