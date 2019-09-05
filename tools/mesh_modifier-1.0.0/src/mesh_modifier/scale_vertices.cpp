
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

namespace mesh_modifier {

void Mesh_modifier::scale_vertices (const std::vector<Point_condition> *pc,
    const Vector<double> center, const Vector<double> scale_vector) {

  // the last unv_container is the one that will be exported
  unsigned uci = unv_container.size() - 1;  

  // find the vertices to be scaled
  std::vector<unsigned int> vertex_index;
  for (int i = 0; i < unv_container[uci].udn_2411.size(); ++i) {
    auto vr = unv_container[uci].udn_2411[i].record2;
    Vector<double> v {vr[0], vr[1], vr[2]};
    bool satisfy_conditions = true;
    for (int j = 0; j < pc->size(); ++j) {
      auto b = pc->at(j).in_condition(v);
      if (!b) {
        satisfy_conditions = false;
        break;
      }
    }
    if (satisfy_conditions) {
      vertex_index.push_back(i);
    }
  }

  // scale them
  for (int i = 0; i < vertex_index.size(); ++i) {
    auto &vr = unv_container[uci].udn_2411[vertex_index[i]].record2;
    Vector<double> v {vr[0], vr[1], vr[2]};
    Vector<double> vsc {(v.x - center.x)*scale_vector.x + center.x,
                        (v.y - center.y)*scale_vector.y + center.y,
                        (v.z - center.z)*scale_vector.z + center.z};
    vr[0] = vsc.x; vr[1] = vsc.y; vr[2] = vsc.z;
  }
} 

} //mesh_modifier
