
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

#include "caviar/objects/shape/polyhedron/output.h"
#include "caviar/objects/shape/polyhedron/polyhedron.h"

#include <string>
#include <cmath>
#include <fstream>

namespace caviar {
namespace objects {
namespace shape {
namespace polyhedron {

Output::Output (CAVIAR *fptr) : Pointers{fptr} {}

Output::~Output () {}

void Output::mesh_povray (const shape::polyhedron::Polyhedron & p_object, std::string file_name) {
  std::ofstream pov_file (file_name.c_str());

    const auto & vertex = p_object.vertex;
    const auto & face = p_object.face;


    pov_file <<"\nmesh {\n";
    for (unsigned int i=0;i<face.size();++i) {
      pov_file << "\ttriangle { ";
      pov_file << "<" << vertex[face[i][0]].x << "," << vertex[face[i][0]].y  << "," << vertex[face[i][0]].z << ">," ;
      pov_file << "<" << vertex[face[i][1]].x << "," << vertex[face[i][1]].y  << "," << vertex[face[i][1]].z << ">,";
      pov_file << "<" << vertex[face[i][2]].x << "," << vertex[face[i][2]].y  << "," << vertex[face[i][2]].z << "> }\n";
    }
    pov_file << "\ttexture {\n\t\tpigment { color rgb<0.9, 0.9, 0.9> }";
    pov_file << "\n\t\tfinish { ambient 0.2 diffuse 0.7 }\n\t}\n";  
    pov_file << "}";
    pov_file << std::flush;

  pov_file.close ();
}

void Output::normals_vectors (const shape::polyhedron::Polyhedron & p_object, std::string file_name) {
  std::ofstream vfptr_file (file_name.c_str());

    const auto & vertex = p_object.vertex;
    const auto & face = p_object.face;
    const auto & normal = p_object.normal;

    //Real_t vec_length = 0.5;
    //Real_t vec_rad = 1.5;

    for (unsigned int i=0;i<face.size();++i) {
      // center of polygons
      Vector<Real_t> cr = { (vertex[face[i][0]].x + vertex[face[i][1]].x + vertex[face[i][2]].x)/3.0, 
                            (vertex[face[i][0]].y + vertex[face[i][1]].y + vertex[face[i][2]].y)/3.0,
                            (vertex[face[i][0]].z + vertex[face[i][1]].z + vertex[face[i][2]].z)/3.0};

      // I think there's two different type of vectors as for output

      // the first type is like this
      // Vector<Real_t> di = cr + normal[i];
      // if one 
      // this Type works for GNUPlot 'splot with vectors'
      Vector<Real_t> di =  normal[i];

      // one can have a scale done on the vectors using gnuplot

      // this line works for when one needs to have the starting point of a
      // normal not to be on the polygon, but on a point a little under it.
      //cr -=  (vec_length/10.0) * normal[i];

      vfptr_file << cr.x << " " << cr.y  << " " << cr.z << " ";
      vfptr_file << di.x << " " << di.y  << " " << di.z << "\n";
      //std::cout << normal[i] << "\n";
    }
    vfptr_file << std::flush;

  vfptr_file.close ();
}

void Output::normals_tcl (const shape::polyhedron::Polyhedron & p_object, std::string file_name) {
  std::ofstream vfptr_file (file_name.c_str());

    const auto & vertex = p_object.vertex;
    const auto & face = p_object.face;
    const auto & normal = p_object.normal;

    Real_t vec_length = 5.0;
    Real_t vec_rad = 1.5;

    for (unsigned int i=0;i<face.size();++i) {
      Vector<Real_t> cr = {(vertex[face[i][0]].x + vertex[face[i][1]].x + vertex[face[i][2]].x)/3.0, 
                            (vertex[face[i][0]].y + vertex[face[i][1]].y + vertex[face[i][2]].y)/3.0,
                            (vertex[face[i][0]].z + vertex[face[i][1]].z + vertex[face[i][2]].z)/3.0};

      Vector<Real_t> di = cr + vec_length*normal[i];
      cr -=  (vec_length/10.0) * normal[i];

      vfptr_file << "graphics top cone ";
      vfptr_file << "{" << cr.x << " " << cr.y  << " " << cr.z << "} ";
      vfptr_file << "{" << di.x << " " << di.y  << " " << di.z << "} ";
      vfptr_file << "radius " << vec_rad << " resolution 5\n";
    }
    vfptr_file << std::flush;

  vfptr_file.close ();
}

void Output::edges_tcl (const shape::polyhedron::Polyhedron & p_object, std::string file_name) {
  std::ofstream vfptr_file (file_name.c_str());

    const auto & vertex = p_object.vertex;
    const auto & edges = p_object.edges;

    Real_t frame_rad = 0.5;
    std::map<std::vector<unsigned int>,std::vector<unsigned int>>::const_iterator it;
    for (it = edges.begin(); it != edges.end(); ++it) {

      vfptr_file << "graphics top cylinder ";
      vfptr_file << "{" << vertex[it->first[0]].x << " " << vertex[it->first[0]].y  << " " << vertex[it->first[0]].z << "} " ;
      vfptr_file << "{" << vertex[it->first[1]].x << " " << vertex[it->first[1]].y  << " " << vertex[it->first[1]].z << "} " ;
      vfptr_file << "radius " << frame_rad << " resolution " << 5 << " filled yes\n";
    }
    vfptr_file << std::flush;

  vfptr_file.close ();
}


void Output::mesh_tcl (const shape::polyhedron::Polyhedron & p_object, std::string file_name) {
  std::ofstream vfptr_file (file_name.c_str());

    const auto & vertex = p_object.vertex;
    const auto & face = p_object.face;

    for (unsigned int i=0;i<face.size();++i) {
      vfptr_file << "graphics top triangle ";
      vfptr_file << "{" << vertex[face[i][0]].x << " " << vertex[face[i][0]].y  << " " << vertex[face[i][0]].z << "} " ;
      vfptr_file << "{" << vertex[face[i][1]].x << " " << vertex[face[i][1]].y  << " " << vertex[face[i][1]].z << "} " ;
      vfptr_file << "{" << vertex[face[i][2]].x << " " << vertex[face[i][2]].y  << " " << vertex[face[i][2]].z << "}\n";
    }

    vfptr_file << std::flush;

  vfptr_file.close ();
}

/*
void Output::mesh_povray (const std::vector<polyhedron::Polyhedron> & shapes) {
  std::ofstream pov_file ("o_mesh.pov");
  for (unsigned int shape_index=0; shape_index < p_object.size();++shape_index) {
    const auto & vertex = shapes[shape_index].vertex;
    const auto & face = shapes[shape_index].face;


    pov_file <<"\nmesh {\n";
    for (unsigned int i=0;i<face.size();++i) {
      pov_file << "\ttriangle { ";
      pov_file << "<" << vertex[face[i][0]].x << "," << vertex[face[i][0]].y  << "," << vertex[face[i][0]].z << ">," ;
      pov_file << "<" << vertex[face[i][1]].x << "," << vertex[face[i][1]].y  << "," << vertex[face[i][1]].z << ">,";
      pov_file << "<" << vertex[face[i][2]].x << "," << vertex[face[i][2]].y  << "," << vertex[face[i][2]].z << "> }\n";
    }
    pov_file << "\ttexture {\n\t\tpigment { color rgb<0.9, 0.9, 0.9> }";
    pov_file << "\n\t\tfinish { ambient 0.2 diffuse 0.7 }\n\t}\n";  
    pov_file << "}";
    pov_file << std::flush;
  }
  pov_file.close ();
}

void Output::normals_vfptr (const std::vector<polyhedron::Polyhedron> & shapes) {
  std::ofstream vfptr_file ("o_normals.tcl");
  for (unsigned int shape_index=0; shape_index < p_object.size();++shape_index) {
    const auto & vertex = shapes[shape_index].vertex;
    const auto & face = shapes[shape_index].face;
    const auto & normal = shapes[shape_index].normal;

    Real_t vec_length = 5.0;
    Real_t vec_rad = 1.5;

    for (unsigned int i=0;i<face.size();++i) {
      Vector<Real_t> cr = {(vertex[face[i][0]].x + vertex[face[i][1]].x + vertex[face[i][2]].x)/3.0, 
                            (vertex[face[i][0]].y + vertex[face[i][1]].y + vertex[face[i][2]].y)/3.0,
                            (vertex[face[i][0]].z + vertex[face[i][1]].z + vertex[face[i][2]].z)/3.0};

      Vector<Real_t> di = cr + vec_length*normal[i];
      cr -=  (vec_length/10.0) * normal[i];

      vfptr_file << "graphics top cone ";
      vfptr_file << "{" << cr.x << " " << cr.y  << " " << cr.z << "} ";
      vfptr_file << "{" << di.x << " " << di.y  << " " << di.z << "} ";
      vfptr_file << "radius " << vec_rad << " resolution 5\n";
    }
    vfptr_file << std::flush;
  }
  vfptr_file.close ();
}

void Output::edges_vfptr (const std::vector<polyhedron::Polyhedron> & shapes) {
  std::ofstream vfptr_file ("o_edges.tcl");
  for (unsigned int shape_index=0; shape_index < p_object.size();++shape_index) {
    const auto & vertex = shapes[shape_index].vertex;
    const auto & edges = shapes[shape_index].edges;

    Real_t frame_rad = 0.5;
    std::map<std::vector<unsigned int>,std::vector<unsigned int>>::const_iterator it;
    for (it = edges.begin(); it != edges.end(); ++it) {

      vfptr_file << "graphics top cylinder ";
      vfptr_file << "{" << vertex[it->first[0]].x << " " << vertex[it->first[0]].y  << " " << vertex[it->first[0]].z << "} " ;
      vfptr_file << "{" << vertex[it->first[1]].x << " " << vertex[it->first[1]].y  << " " << vertex[it->first[1]].z << "} " ;
      vfptr_file << "radius " << frame_rad << " resolution " << 5 << " filled yes\n";
    }
    vfptr_file << std::flush;
  }
  vfptr_file.close ();
}


void Output::mesh_vfptr (const std::vector<polyhedron::Polyhedron> & shapes) {
  std::ofstream vfptr_file ("o_mesh.tcl");
  for (unsigned int shape_index=0; shape_index < p_object.size();++shape_index) {
    const auto & vertex = shapes[shape_index].vertex;
    const auto & face = shapes[shape_index].face;

    for (unsigned int i=0;i<face.size();++i) {
      vfptr_file << "graphics top triangle ";
      vfptr_file << "{" << vertex[face[i][0]].x << " " << vertex[face[i][0]].y  << " " << vertex[face[i][0]].z << "} " ;
      vfptr_file << "{" << vertex[face[i][1]].x << " " << vertex[face[i][1]].y  << " " << vertex[face[i][1]].z << "} " ;
      vfptr_file << "{" << vertex[face[i][2]].x << " " << vertex[face[i][2]].y  << " " << vertex[face[i][2]].z << "}\n";
    }

    vfptr_file << std::flush;
  }
  vfptr_file.close ();
}

*/
} //polyhedron
} //shape
} //objects
} // namespace caviar

