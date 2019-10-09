
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

#include "caviar/objects/shape/polyhedron/preprocess.h"
#include "caviar/objects/shape/polyhedron/polyhedron.h"
#include "caviar/interpreter/error.h"

#include <string>
#include <cmath>
#include <fstream>

namespace caviar {
namespace objects {
namespace shape {
namespace polyhedron {

Preprocess::Preprocess (CAVIAR *fptr) : Pointers{fptr}
 {}

Preprocess::~Preprocess () { }

// This function makes the 'vertex_map'
void Preprocess::merge_vertices (shape::polyhedron::Polyhedron & p_object) { // There can be a function that merge vertices which are closer than a small distance.
  auto & vertex = p_object.vertex;
  auto & vertex_map = p_object.vertex_map;

  vertex_map.resize (vertex.size());
  for (unsigned int i=0;i<vertex.size();++i) {
    for (unsigned int j=i+1;j<vertex.size();++j) {
       if (vertex[i]==vertex[j]) {
        vertex_map[j].push_back(i);
      }
    }
  }
}

// It's hard to explane what this function does part by part. What matters is 
// that this function changes the order of the triangles vertices so that when
// the normal vectors are made, all of them be in the same direction.
void Preprocess::pre_correct_normals (shape::polyhedron::Polyhedron & p_object) {

  auto & face = p_object.face;
  const auto & edges = p_object.edges;

  if (face.size() == 0 || edges.size() == 0) error-> all(FC_FILE_LINE_FUNC,"The polyhedron face or edge is empty");

  std::vector<unsigned int> face_list; // sequence of faces
  std::vector<unsigned int> face_in,face_out; // face index of edge_list
  std::vector<std::vector<unsigned int>> edge_list; // sequence of edges

  int start_face = 0;

  face_list.push_back (start_face);
  std::cout << "edges.size: "<< edges.size() << "\n";
  std::cout << "face_list.size: "<< face_list.size() << "\n";
  std::cout << "face.size: "<< face.size() << "\n";  
  std::cout << "edge_list.size: "<< edge_list.size() << "\n";
std::cout << "XXX HEY 1 " << std::endl;
//   while (face_list.size() != face.size()) 
  {
    for (unsigned int i = 0; i<face_list.size(); ++i) {

      std::vector<unsigned int> e1,e2,e3;

#define MAKE_EDGE(NUM1,NUM2,EN1)                                     \
      {                                                              \
        int v1 = face[face_list[i]][NUM1];                           \
        int v2 = face[face_list[i]][NUM2];                           \
        if (v1<v2) {EN1.push_back(v1); EN1.push_back(v2);}           \
        else {EN1.push_back(v2); EN1.push_back(v1);}                 \
      }

      MAKE_EDGE(0,1,e1);

      MAKE_EDGE(1,2,e2);

      MAKE_EDGE(2,0,e3);

#undef MAKE_EDGE

      std::map<std::vector<unsigned int>,std::vector<unsigned int>>::const_iterator it;    
#define MAKE_EDGE_LIST(EN1)                                        \
      {                                                            \
        it = edges.find (EN1);                                     \
        if (it == edges.end() )                                    \
          error->all (FC_FILE_LINE_FUNC, "a bug in the code 3!");  \
        if (it->second.size() == 2) { \
        unsigned int its01;                                        \
        if (it->second[0] == face_list[i]) its01 = it->second[1];  \
        else its01 = it->second[0];                                \
        bool face_used = false;                                    \
        for (unsigned int j=0;j<face_list.size();++j) {            \
          if (its01 == face_list[j]) {                             \
            face_used = true; break;                               \
          }                                                        \
        }                                                          \
        if (!face_used) {                                          \
          face_list.push_back (its01);                             \
          edge_list.push_back (EN1);                               \
          face_in.push_back (face_list[i]);                        \
          face_out.push_back (its01);                              \
        }                                                          \
        } \
      }
      MAKE_EDGE_LIST(e1);
      MAKE_EDGE_LIST(e2);
      MAKE_EDGE_LIST(e3);
#undef MAKE_EDGE_LIST

    }
  }

  //std::cout << "edges.size: "<< edges.size() << "\n";
  //std::cout << "face_list.size: "<< face_list.size() << "\n";
  //std::cout << "face.size: "<< face.size() << "\n";  
  //std::cout << "edge_list.size: "<< edge_list.size() << "\n";
std::cout << "XXX HEY 2 " << std::endl;
  for (unsigned int i=0; i<edge_list.size(); ++i) {

      std::map<std::vector<unsigned int>,std::vector<unsigned int>>::const_iterator it;    
      it = edges.find (edge_list[i]);  
      if (it == edges.end())   {error->all (FC_FILE_LINE_FUNC, "a bug in the code 1!");  }//
      if (it->second.size() != 2){ error->all (FC_FILE_LINE_FUNC, "a bug in the code 2!");  }//

      {
          unsigned int v1 = it->first[0]; // vertex index
          unsigned int v2 = it->first[1]; // 
          unsigned int f1 = it->second[0];// face index
          unsigned int f2 = it->second[1];

          bool f1_right,f2_right ; // true or false whether face1 or face2 edges are ordered in a right hand twist order.
          int f1_case=-1,f2_case=-1; // there's three different right hand and three left hand cases in which the vertices are ordered.
          if (v1==face[f1][0]) {
            if (v2==face[f1][1]) {
              f1_right = true;
              f1_case = 1;
            } else { // (v2==face[f1][2]) 
              f1_right = false;
              f1_case = 2;
            }
          } else if (v1==face[f1][1]) {
            if (v2==face[f1][0]) {
              f1_right = false;
              f1_case = 3;
            } else { // (v2==face[f1][2]) 
              f1_right = true;
              f1_case = 4;
            }
          } else { // (v1==face[f1][2])
            if (v2==face[f1][0]) {
              f1_right = true;
              f1_case = 5;
            } else { // (v2==face[f1][1]) 
              f1_right = false;
              f1_case = 6;
            } 
          }

          if (v1==face[f2][0]) {
            if (v2==face[f2][1]) {
              f2_right = true;
              f2_case = 1;
            } else { // (v2==face[f2][2]) 
              f2_right = false;
              f2_case = 2;
            }
          } else if (v1==face[f2][1]) {
              if (v2==face[f2][0]) {
                f2_right = false;
                f2_case = 3;
            } else { // (v2==face[f2][2]) 
              f2_right = true;
              f2_case = 4;
            }
          } else { // (v1==face[f2][2])
            if (v2==face[f2][0]) {
              f2_right = true;
              f2_case = 5;
            } else { // (v2==face[f2][1]) 
              f2_right = false;
              f2_case = 6;
            } 
          }

          if (f1_right==f2_right) { // if two neighborlist faces has similar twist order in one edge, 
                                    // their normals will be directed differently. so it has to be corrected.
            if (f2 == face_out[i]) { // only change the outer face

              switch (f2_case) {
                case 1: case 3: {
                  unsigned int tmp = face[f2][0];
                  face[f2][0] = face[f2][1];
                  face[f2][1] = tmp; 
                }  break;
                case 2: case 5: {
                  unsigned int tmp = face[f2][0];
                  face[f2][0] = face[f2][2];
                  face[f2][2] = tmp;
                } break;
                case 4: case 6: {
                  unsigned int tmp = face[f2][1];
                  face[f2][1] = face[f2][2];
                  face[f2][2] = tmp;
                }  break;
              }
            } else {

              switch (f1_case) {
                case 1: case 3: {
                  unsigned int tmp = face[f1][0];
                  face[f1][0] = face[f1][1];
                  face[f1][1] = tmp; 
                } break;
                case 2: case 5: {
                  unsigned int tmp = face[f1][0];
                  face[f1][0] = face[f1][2];
                  face[f1][2] = tmp;
                } break;
                case 4: case 6: {
                  unsigned int tmp = face[f1][1];
                  face[f1][1] = face[f1][2];
                  face[f1][2] = tmp;
                }  break;
              }
            }

        }
      }
  }
}


} //polyhedron
} //shape
} //objects
} // namespace caviar

