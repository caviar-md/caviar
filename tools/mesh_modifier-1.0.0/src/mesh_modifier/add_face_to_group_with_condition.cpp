
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

void Mesh_modifier::add_face_to_group_with_condition (const std::string gr, const std::vector<Point_condition> *pc) {

  // the last unv_container is the one that will be exported
  unsigned uci = unv_container.size() - 1;  

  // find and set the group_index if a group with the name 'gr' exists
  int group_index = 0;
  bool new_group = true;
  for (int i = 0; i<unv_container[uci].udn_2467.size(); ++i) {
    auto record2 =  unv_container[uci].udn_2467[i].record2;
   if (std::to_string(record2)==gr) {
      new_group = false;
      group_index = i;
      break;
   }
  }

  // create a new group if 'gr' is not found.
  if (new_group) {
    Universal_dataset_number_2467 u;
    u.record1.resize(8);
    u.record1[0] = unv_container[uci].udn_2467.size() + 1;
    u.record1[1] = 0; u.record1[2] = 0; u.record1[3] = 0; u.record1[4] = 0;
    u.record1[5] = 0; u.record1[6] = 0; u.record1[7] = 0;
    u.record2 = std::stoi(gr); // XXX for now 'record2' is of int type
    unv_container[uci].udn_2467.push_back(u);
    group_index = unv_container[uci].udn_2467.size() - 1;
  } 

  auto &u2467 = unv_container[uci].udn_2467[group_index];

  auto one_third = 1/3.0;

  // find, check and add the face elements to the groups
  for (int i = 0; i < unv_container[uci].udn_2412.size(); ++i) {

    auto &u2412 = unv_container[uci].udn_2412[i];

    // here we find face center
    Vector<double> face_center;
    // plane triangle cases
    if (u2412.record1[1] == 41) {
      auto v1 = unv_container[uci].udn_2411[u2412.record2[0]-1].record2;
      auto v2 = unv_container[uci].udn_2411[u2412.record2[1]-1].record2;
      auto v3 = unv_container[uci].udn_2411[u2412.record2[2]-1].record2;
      Vector<double> vv1 {v1[0], v1[1], v1[2]};
      Vector<double> vv2 {v2[0], v2[1], v2[2]};
      Vector<double> vv3 {v3[0], v3[1], v3[2]};
      face_center = one_third * (vv1 + vv2 + vv3);

    // plane Linear Quadrilateral
    } else if (u2412.record1[1] == 44) {

      auto v1 = unv_container[uci].udn_2411[u2412.record2[0]-1].record2;
      auto v2 = unv_container[uci].udn_2411[u2412.record2[1]-1].record2;
      auto v3 = unv_container[uci].udn_2411[u2412.record2[2]-1].record2;
      auto v4 = unv_container[uci].udn_2411[u2412.record2[3]-1].record2;
      Vector<double> vv1 {v1[0], v1[1], v1[2]};
      Vector<double> vv2 {v2[0], v2[1], v2[2]};
      Vector<double> vv3 {v3[0], v3[1], v3[2]};
      Vector<double> vv4 {v4[0], v4[1], v4[2]};
      face_center = 0.25 * (vv1 + vv2 + vv3 + vv4);

    } else {
      // do nothing yet.
      continue;
    }


    // now we check the face_center whether it satisfy all the conditions:
    bool satisfy_conditions = true;
    for (int j = 0; j < pc->size(); ++j) {
      auto b = pc->at(j).in_condition(face_center);
      if (!b) {
        satisfy_conditions = false;
        break;
      }
    }

    // now add the face element to the group
    if (satisfy_conditions) {
      int element_label = u2412.record1[0];
      u2467.record3.push_back(8);
      u2467.record3.push_back(element_label); 
      u2467.record3.push_back(0);
      u2467.record3.push_back(0);
    }

  }

  // setting or re-setting the number of elements in the group.
  u2467.record1[7] = static_cast<int> (u2467.record3.size() / 4);

} 

} //mesh_modifier
