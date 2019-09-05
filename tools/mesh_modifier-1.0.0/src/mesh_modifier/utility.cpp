
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

#include "mesh_modifier.h"

#include <iostream>
#include <algorithm>    // std::sort, std::count

namespace mesh_modifier {

  void Mesh_modifier::remove_edge (std::vector<unsigned> &edge_k
                                         ,std::vector<unsigned> &removed_entity_label) {

    unsigned uci = unv_container.size() - 1;
    auto & udn_2412 = unv_container[uci].udn_2412;
    std::vector<unsigned> to_remove;
    for (unsigned i = 0 ; i < udn_2412.size(); ++i) {
      if (udn_2412[i].record1[1]==11) {
        std::vector<unsigned> e {udn_2412[i].record3[0], udn_2412[i].record3[1]};
        if (e[0] > e[1]) {
          unsigned tmp=e[0]; e[0]=e[1]; e[1]=tmp;
        }
        
        if (e == edge_k) {
          to_remove.push_back(i);
          removed_entity_label.push_back(udn_2412[i].record1[0]);          
        }
      }
    }
    std::sort(to_remove.begin(), to_remove.end(), std::greater<unsigned>());
    for (unsigned i = 0 ; i < to_remove.size(); ++i) {
      udn_2412.erase(udn_2412.begin() + to_remove[i]);
    }
  }
  
  
  void Mesh_modifier::remove_quad (std::vector<unsigned> &quad_k
                                         ,std::vector<unsigned> &removed_entity_label) {

    unsigned uci = unv_container.size() - 1;
    auto & udn_2412 = unv_container[uci].udn_2412;
    std::vector<unsigned> to_remove;
    for (unsigned i = 0 ; i < udn_2412.size(); ++i) {
      if (udn_2412[i].record1[1]==44) {
        std::vector<unsigned> e {udn_2412[i].record2[0], udn_2412[i].record2[1]
                                ,udn_2412[i].record2[2], udn_2412[i].record2[3]};
                 
        std::sort(e.begin(), e.end());
        
        if (e == quad_k) {
          to_remove.push_back(i);
          removed_entity_label.push_back(udn_2412[i].record1[0]);
        }
      }
    }

    std::sort(to_remove.begin(), to_remove.end(), std::greater<unsigned>());
    for (unsigned i = 0 ; i < to_remove.size(); ++i) {
      udn_2412.erase(udn_2412.begin() + to_remove[i]);
      
    }  

  }
  
  void Mesh_modifier::add_edge_to_vector (std::vector<std::vector<unsigned>> &hexa_edges
                                  ,unsigned pl1, unsigned pl2) {
    std::vector<unsigned> e {pl1, pl2};
    std::sort (e.begin(), e.end());
    hexa_edges.push_back (e);
  }

  void Mesh_modifier::add_triangle_to_vector (std::vector<std::vector<unsigned>> &tetra_triangle
                                  ,unsigned pl1, unsigned pl2
                                  ,unsigned pl3) {
    std::vector<unsigned> e {pl1, pl2, pl3};
    std::sort (e.begin(), e.end());
    tetra_triangle.push_back (e);
  }
  
  void Mesh_modifier::add_quad_to_vector (std::vector<std::vector<unsigned>> &hexa_quads
                                  ,unsigned pl1, unsigned pl2
                                  ,unsigned pl3, unsigned pl4) {
    std::vector<unsigned> e {pl1, pl2, pl3, pl4};
    std::sort (e.begin(), e.end());
    hexa_quads.push_back (e);
  }
  
    void Mesh_modifier::add_tetrahedron_to_vector (std::vector<std::vector<unsigned>> &tetrahedrons
                                  ,unsigned pl1, unsigned pl2
                                  ,unsigned pl3, unsigned pl4) {
    std::vector<unsigned> e {pl1, pl2, pl3, pl4};
    std::sort (e.begin(), e.end());
    tetrahedrons.push_back (e);
  }
}
