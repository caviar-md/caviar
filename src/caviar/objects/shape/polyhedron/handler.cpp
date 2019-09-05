
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

#include "caviar/objects/shape/polyhedron/handler.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/shape/polyhedron/polyhedron.h"
#include "caviar/objects/shape/polyhedron/input.h"
#include "caviar/objects/shape/polyhedron/utility.h"
#include "caviar/objects/shape/polyhedron/preprocess.h"
#include "caviar/objects/shape/polyhedron/postprocess.h"
#include "caviar/objects/shape/polyhedron/point_inside.h"
#include "caviar/objects/shape/polyhedron/output.h"

#include "caviar/objects/shape/polyhedron/format_vtk_reader.h"

#include <string>
#include <cmath>
#include <fstream>

namespace caviar {
namespace objects {
namespace shape {
namespace polyhedron {

Handler::Handler (CAVIAR *fptr) : Pointers{fptr},
  polyhedron_input {new shape::polyhedron::Input{fptr}},
  polyhedron_preprocess {new shape::polyhedron::Preprocess{fptr}},
  polyhedron_postprocess {new shape::polyhedron::Postprocess{fptr}},
  polyhedron_utility {new shape::polyhedron::Utility{fptr}},
  polyhedron_point_inside {new shape::polyhedron::Point_Inside {fptr}},
  polyhedron_output {new shape::polyhedron::Output {fptr}},  
  polyhedron_read{false},
  output_mesh_tcl{false},
  output_normals_tcl{false},
  output_edges_tcl{false},
  output_mesh_povray{false},
  output_normals_vectors{false},
  invert_normals{false},
  correct_normals{false},
  use_grid{false},
  an_inside_point_is_set{false}

  {
    radius_max = -1.0;
  }



Handler::~Handler () {
  delete polyhedron_input;
  delete polyhedron_preprocess;
  delete polyhedron_postprocess;
  delete polyhedron_utility;
  delete polyhedron_point_inside;
  delete polyhedron_output;  
}

bool Handler::read(class caviar::interpreter::Parser * parser) {
  output->info("Polyhedron read:");
  bool in_file = true;
      
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    const auto t = token.string_value;
    if (string_cmp(t,"vtk_file_name")) {
      const auto token = parser->get_val_token();
      const auto file_name = token.string_value;
      polyhedron_input -> read_vtk (polyhedron, file_name);
      polyhedron_read = true;
    } else if (string_cmp(t,"stl_file_name")) {
      const auto token = parser->get_val_token();
      const auto file_name = token.string_value;
      polyhedron_input -> read_stl (polyhedron, file_name);
      polyhedron_read = true;
    } else if (string_cmp(t,"output")) {
      const auto token = parser->get_val_token();
      const auto t = token.string_value; 
      if (string_cmp(t,"vectors") || string_cmp(t,"vector")) {
        output_normals_vectors = true;
      }

      if (string_cmp(t,"vfptr") || string_cmp(t,"tcl")) {
        output_mesh_tcl = true;
        output_normals_tcl = true;
        output_edges_tcl = true;
      }
    } else if (string_cmp(t,"generate")) {
      command_generate ();
    } else if (string_cmp(t,"use_grid")) {
      use_grid = true; return true;
    } else if (string_cmp(t,"an_inside_point")) {
      GET_OR_CHOOSE_A_REAL(an_inside_point.x,"","")
      GET_OR_CHOOSE_A_REAL(an_inside_point.y,"","")
      GET_OR_CHOOSE_A_REAL(an_inside_point.z,"","")
      an_inside_point_is_set = true;      
    } else if (string_cmp(t,"write_unstructured_vtk")) {
      class Format_vtk_reader fvr (fptr);
      fvr.write_unstructured_vtk4 (polyhedron);      
    } else if (string_cmp(t,"write_polydata_vtk")) {
      class Format_vtk_reader fvr (fptr);
      fvr.write_polydata_vtk4 (polyhedron);      
    } else if (string_cmp(t,"point_is_inside_method")) {
      int i = 0;
      GET_OR_CHOOSE_A_INT(i,"","")
      polyhedron_point_inside->point_is_inside_method = i; 
    } else if (string_cmp(t,"thickness")) {
      auto thickness_ = parser->get_literal_real();
      polyhedron.thickness = thickness_;
    } else if (string_cmp(t,"radius_max")) {
      radius_max = parser->get_literal_real();
      polyhedron.grid_tol = radius_max;
    } else if (string_cmp(t,"invert_normals")) {
      invert_normals = true;
    } else if (string_cmp(t,"correct_normals")) {
      correct_normals = true;
    } else if (string_cmp(t,"grid")) {
      auto nx_ = parser->get_literal_int();
      auto ny_ = parser->get_literal_int();
      auto nz_ = parser->get_literal_int();
      if (nx_ < 1 || ny_ < 1 || nz_ < 1) error->all(FC_FILE_LINE_FUNC_PARSE, "grids has to be larger than 1");
      polyhedron.nx_part = nx_;     polyhedron.nx_part = ny_;    polyhedron.nx_part = nz_;
    } else error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
  }
  
  return in_file;;
}

bool Handler::is_inside (const Vector<double> &v0) {
  return  polyhedron_point_inside -> is_inside (polyhedron, v0);   
}


bool Handler::is_inside (const Vector<double> &v, const double r) {
  if (use_grid)
    return polyhedron_point_inside -> is_inside_grid (polyhedron, v, r); 

  return polyhedron_point_inside -> is_inside_all (polyhedron, v, r); 
}

bool Handler::in_contact (const Vector<double> &v, const double r, Vector<double> & contact_vector) {
  if (use_grid)
    return polyhedron_point_inside -> in_contact_grid (polyhedron, v, r, contact_vector);

  return polyhedron_point_inside -> in_contact_all    (polyhedron, v, r, contact_vector);
}



void  Handler::command_generate () {
  if (!polyhedron_read) {    
      error->all (FC_FILE_LINE_FUNC, "No polyhedron is imported.");      
  }

  if (correct_normals) {    
    output->info("polyhedron: generate : correct_normals");
    polyhedron_preprocess -> pre_correct_normals (polyhedron);
  }
    
  output->info("polyhedron: generate : make_normal");
  polyhedron_utility -> make_normal (polyhedron);

  output->info("polyhedron: generate : make_edge_norms");
  polyhedron_utility -> make_edge_norms (polyhedron);
  
  if (an_inside_point_is_set) {
    invert_normals = polyhedron_utility -> normals_are_pointing_outside(polyhedron, an_inside_point);
  }

  if (invert_normals) {
    output->info("polyhedron: generate : invert_normals");
    polyhedron_utility -> invert_normals (polyhedron);
  }
  
  output->info("polyhedron: generate : lowest_highest_coord : ",0 , false);
  polyhedron_postprocess -> lowest_highest_coord (polyhedron);

  if (use_grid) {
    if (radius_max < 0) 
      error->all (FC_FILE_LINE_FUNC, "Radius_max is not set.");      

    output->info("polyhedron: generate : make_grid");
    polyhedron_postprocess -> make_grid (polyhedron);
  }    


  if (output_mesh_povray)  polyhedron_output -> mesh_povray (polyhedron); 
  if (output_mesh_tcl)  polyhedron_output -> mesh_tcl (polyhedron); 
  if (output_normals_tcl)  polyhedron_output -> normals_tcl (polyhedron);
  if (output_normals_vectors)  polyhedron_output -> normals_vectors (polyhedron);
  if (output_edges_tcl)  polyhedron_output -> edges_tcl (polyhedron);

}

} //polyhedron
} //shape
} //objects
} // namespace caviar

