
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

#include "caviar/objects/shape/polyhedron/utility.h"
#include "caviar/objects/shape/polyhedron/polyhedron.h"
#include "caviar/utility/interpreter_io_headers.h"

#include <string>
#include <cmath>
#include <fstream>


CAVIAR_NAMESPACE_OPEN

namespace shape {
namespace polyhedron {

Utility::Utility (CAVIAR *fptr) : Pointers{fptr} {}

Utility::~Utility () { }


void Utility::invert_normals (shape::polyhedron::Polyhedron & p_object) {
  auto & normal = p_object.normal;
  for (unsigned int i=0;i<normal.size();++i)
    normal[i] *= -1;
}


void Utility::make_edge_norms (shape::polyhedron::Polyhedron & p_object) {
  const auto & vertex = p_object.vertex;
  const auto & face = p_object.face;
  const auto & edges = p_object.edges;
  const auto & normal = p_object.normal;
  auto & edge_norms3 = p_object.edge_norms3;

  std::vector<std::vector<Vector<Real_t>>> edge_norms1, edge_norms2;
  std::map<std::vector<unsigned int>,std::vector<unsigned int>>::const_iterator it_edges; 
  const unsigned int fsize = face.size();
  edge_norms1.resize (fsize);  
  edge_norms2.resize (fsize);  
  edge_norms3.resize (fsize);
  for (unsigned int i=0;i<fsize;++i) {
    for (unsigned int j=0;j<face[i].size();++j) {
      Vector<Real_t> n1 = normal[i];
      unsigned int k0 = j; unsigned int k1 = (j==face[i].size()-1) ? 0 : j+1;
      unsigned int v0 = face[i][k0]; unsigned int v1 = face[i][k1];
      if (v0<v1) 
        it_edges = edges.find (std::vector<unsigned int>{v0,v1});
      else
        it_edges = edges.find (std::vector<unsigned int>{v1,v0});
      if (it_edges->second.size() == 2) {
        if (it_edges->second[0] != i )
          n1 += normal [it_edges->second[0]];
        else
          n1 += normal [it_edges->second[1]];
      }
      n1 /= std::sqrt (n1*n1); // it is not necessary
      edge_norms1[i].push_back (n1);
//-------
      Vector<Real_t> n2 = vertex[face[i][k1]] - vertex[face[i][k0]];    
      n2 /= std::sqrt (n2*n2); // it is not necessary
      edge_norms2[i].push_back (n2);
//-------
      Vector<Real_t> n3 = cross_product (n1,n2);
      n3 /= std::sqrt (n3*n3); // it is not necessary
      edge_norms3[i].push_back (n3);    
    }
  }


}



void Utility::make_normal (shape::polyhedron::Polyhedron & p_object) { // It's supposed here that the vertices are written in order (right-hand rotation)
  const auto & vertex = p_object.vertex;
  const auto & face = p_object.face;
  auto & normal = p_object.normal;

  for (unsigned int i=0;i<face.size();++i) {
    if (face[i].size()<3) continue;
    //std::cout << i << " , v1: 1 " << face[i][1] << std::endl;
    //std::cout << i << " , v1: 2 " << face[i][2] << std::endl;
    //std::cout << i << " , v1: 0 " << face[i][0] << std::endl;
    Vector<Real_t> v1 = vertex[face[i][1]] - vertex[face[i][0]];
    Vector<Real_t> v2 = vertex[face[i][2]] - vertex[face[i][0]];
    Vector<Real_t> n = cross_product (v1,v2);
//    std::cout<<"v1: "<<v1 << " v2 :" << v2 << " v1*v2: "<< n <<std::endl;
    Real_t n_lenght = sqrt (n.x*n.x + n.y*n.y + n.z*n.z);
    n /= n_lenght;
    normal.push_back(n);
  }

}

bool Utility::normals_are_pointing_outside(shape::polyhedron::Polyhedron & p_object, const Vector<double> &v) {
  const auto & vertex = p_object.vertex;
  const auto & face = p_object.face;
  auto & normal = p_object.normal;

  double nearest_distance=1e9;
  int nearest_index = -1;
  // find the nearest face center;
  for (unsigned int i=0;i<face.size();++i) {

    // center of the faces
    Vector<Real_t> cr = { (vertex[face[i][0]].x + vertex[face[i][1]].x + vertex[face[i][2]].x)/3.0, 
                          (vertex[face[i][0]].y + vertex[face[i][1]].y + vertex[face[i][2]].y)/3.0,
                          (vertex[face[i][0]].z + vertex[face[i][1]].z + vertex[face[i][2]].z)/3.0};

    auto dis_sq = (v-cr)*(v-cr);
    
    if (nearest_distance > dis_sq) {
      nearest_distance = dis_sq;
      nearest_index = i;
    }
  }


  if (nearest_index==-1)
    error->all(FC_FILE_LINE_FUNC,"'An_inside_point' is far from any polygons on the geometry");

  unsigned int i = nearest_index;
  Vector<Real_t> cr = { (vertex[face[i][0]].x + vertex[face[i][1]].x + vertex[face[i][2]].x)/3.0, 
                        (vertex[face[i][0]].y + vertex[face[i][1]].y + vertex[face[i][2]].y)/3.0,
                        (vertex[face[i][0]].z + vertex[face[i][1]].z + vertex[face[i][2]].z)/3.0};

  double dot_product = normal[i] * (v-cr) ;
  if ( dot_product == 0.0 )
    error->all(FC_FILE_LINE_FUNC,"'An_inside_point' gives a zero dot_product. Change the point");

  if (dot_product < 0.0) 
    return true;
  else
    return false;  
}

} //polyhedron
} //shape

CAVIAR_NAMESPACE_CLOSE


