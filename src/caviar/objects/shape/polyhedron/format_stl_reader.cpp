
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

#include "caviar/objects/shape/polyhedron/format_stl_reader.h"
#include "caviar/objects/shape/polyhedron/polyhedron.h"
#include "caviar/objects/shape/polyhedron/preprocess.h"
#include "caviar/utility/interpreter_io_headers.h"

namespace caviar {

namespace shape {
namespace polyhedron {

Format_stl_reader::Format_stl_reader (CAVIAR *fptr) : Pointers{fptr} {}

Format_stl_reader::~Format_stl_reader() {}

void Format_stl_reader::read_polyhedron (shape::polyhedron::Polyhedron &p_object, const std::string & file) {
  class caviar::interpreter::Parser *parser = new caviar::interpreter::Parser {fptr, file};
  std::cout << "info: Stl_file: reading a Stl_file Stl file: "<< file << std::endl;

  auto & vertex = p_object.vertex;
  auto & face = p_object.face;
  auto & edges = p_object.edges;
  auto & face_id = p_object.face_id;

  if (face.size() != face_id.size()) {
    error->all(FC_FILE_LINE_FUNC,"face.size() != face_id.size()");
  }

  //auto & vertex_map = p_object.vertex_map;
  std::cout << "info: Stl_file: recording Stl file as a 3D shape with index " << " ? " << std::endl;


  bool endsolid_found = false;
  int no_polygons = 0;
  while(true) {
    auto t = parser->get_val_token ();
    auto ts = t.string_value;

    if (t.kind==caviar::interpreter::Kind::eof) break;

    if (t.kind==caviar::interpreter::Kind::identifier) {
      if (ts == "solid") {
        parser->go_to_next_line();
        continue;
      }

      if (ts == "endsolid") {
        endsolid_found = true;
        break;
      }

      if (ts == "facet") {
        ++no_polygons;

        t = parser->get_val_token (); 
        ts = t.string_value;
        if (ts != "normal") 
          error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown stl format: expected 'normal' string");

        double nx = parser->get_literal_real (); 
        double ny = parser->get_literal_real (); 
        double nz = parser->get_literal_real (); 

        Vector<double> normal {nx, ny, nz}; //XXX
        normal.x = normal.x;

        parser->go_to_next_line();
        continue;        
        
      }

      if (ts == "outer") {
        // initial string token
        t = parser->get_val_token (); 
        ts = t.string_value;
        if (ts != "loop") 
          error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown stl format: expected 'loop' string");

        parser->go_to_next_line();

        // getting the vertex data from the 'stl' file
        t = parser->get_val_token (); 
        ts = t.string_value;
        if (ts != "vertex") 
          error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown stl format: expected 'vertex' string");

        double v1x = parser->get_literal_real (); 
        double v1y = parser->get_literal_real (); 
        double v1z = parser->get_literal_real (); 

        Vector<double> v1 {v1x, v1y, v1z};

        parser->go_to_next_line();

        t = parser->get_val_token (); 
        ts = t.string_value;
        if (ts != "vertex") 
          error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown stl format: expected 'vertex' string");

        double v2x = parser->get_literal_real (); 
        double v2y = parser->get_literal_real (); 
        double v2z = parser->get_literal_real (); 

        Vector<double> v2 {v2x, v2y, v2z};

        parser->go_to_next_line();

        t = parser->get_val_token (); 
        ts = t.string_value;
        if (ts != "vertex") 
          error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown stl format: expected 'vertex' string");

        double v3x = parser->get_literal_real (); 
        double v3y = parser->get_literal_real (); 
        double v3z = parser->get_literal_real (); 

        Vector<double> v3 {v3x, v3y, v3z};

        parser->go_to_next_line();

        // setting the data into the Polyhedron instance

        vertex.push_back (v1);
        vertex.push_back (v2);
        vertex.push_back (v3);


        // final string tokens
        t = parser->get_val_token (); 
        ts = t.string_value;
        if (ts != "endloop") 
          error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown stl format: expected 'endloop' string");

        parser->go_to_next_line();

        t = parser->get_val_token (); 
        ts = t.string_value;
        if (ts != "endfacet") 
          error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown stl format: expected 'facet' string");

        parser->go_to_next_line();
        continue;
      }
    }

  }

  if (!endsolid_found) {
    output->warning("missing 'endsolid'. Maybe an incomplete stl file.");
  }

  // merge vertices and make the vertex_map out of it.
  polyhedron::Preprocess p_pre (fptr);
  p_pre.merge_vertices(p_object); 

  for (unsigned int i=0;i<static_cast<unsigned int>(no_polygons);++i) {

    unsigned int num1 = 3*i;
    unsigned int num2 = num1 + 1;
    unsigned int num3 = num1 + 2;


    //num1 = (vertex_map[num1].size()==0)?num1 :vertex_map[num1][1];  // uses vertex_map instead of vertex
    //num2 = (vertex_map[num2].size()==0)?num2 :vertex_map[num2][1];  //
    //num3 = (vertex_map[num3].size()==0)?num3 :vertex_map[num3][1];  //

//    std::cout << "f: " << num1 << " , " << num2 << " , " << num3 << std::endl;

    std::vector<unsigned int> gons;
    gons.push_back (num1); gons.push_back (num2); gons.push_back (num3);
    face.push_back (gons);      
    face_id.push_back(-1); // invalid face_id;

    //faces_of_vertex[num1].push_back(i); faces_of_vertex[num2].push_back(i); faces_of_vertex[num3].push_back(i);

    std::map<std::vector<unsigned int>,std::vector<unsigned int>>::iterator it_edges; 
#define ADD_EDGE(NUM1,NUM2)                                    \
    {                                                      \
      std::vector<unsigned int> check_face = {i};          \
      std::vector<unsigned int> temp_edge;                 \
      if (NUM1<NUM2) temp_edge = {NUM1,NUM2};              \
      else temp_edge = {NUM2,NUM1};                        \
      it_edges = edges.find(temp_edge);                    \
      if (it_edges == edges.end()) {                       \
        edges.insert (make_pair(temp_edge,check_face));    \
      }  else {                                            \
        it_edges->second.push_back(i);                     \
      }                                                    \
    }

    ADD_EDGE(num1,num2);
    ADD_EDGE(num2,num3);
    ADD_EDGE(num3,num1);

#undef ADD_EDGE

  }
      
  delete parser;
}

} //polyhedron
} //shape

} // namespace caviar

