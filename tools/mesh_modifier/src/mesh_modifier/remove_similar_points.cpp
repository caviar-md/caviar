
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
#include <algorithm> // std::sort, std::count

#define POW2(X) X *X
namespace mesh_modifier
{

  // the usage of similar points in a geometry is when one wants to move an element
  // even when it has some shared vertex position with other elements.
  void Mesh_modifier::remove_similar_points(const double tol_sqr)
  {
    unsigned uci = unv_container.size() - 1;

    auto &udn_2411 = unv_container[uci].udn_2411;
    auto num_of_points = udn_2411.size();

    // this part is used to refix the point id sequence from 1 to the last number
    // after some of them are deleted.
    // this may looks unnecessary for unv format, but since we need the record in
    // the exporting data process in other formats such as vtk and stl, fixing it
    // here prevents unexpected problems and bugs.
    std::vector<int> new_record1_0;

    // there can be more than one similar point to one position, so we make a
    // list of them
    std::vector<std::vector<unsigned>> similar_points;

    std::vector<unsigned> similar_list;
    std::vector<unsigned> to_remove;
    similar_points.resize(num_of_points);

    new_record1_0.resize(num_of_points + 1, 0);
    for (unsigned i = 0; i < num_of_points; ++i)
    {

      int id = udn_2411[i].record1[0];
      new_record1_0[id] = id;

      // if number of (i+1) elements in 'similar_list' is positive, it means that
      // it has been used as a 'similar point' to be deleted.
      if (std::count(similar_list.begin(), similar_list.end(), id) > 0)
        continue;

      //
      for (unsigned j = i + 1; j < num_of_points; ++j)
      {
        double distance = POW2((udn_2411[i].record2[0] - udn_2411[j].record2[0])) + POW2((udn_2411[i].record2[1] - udn_2411[j].record2[1])) + POW2((udn_2411[i].record2[2] - udn_2411[j].record2[2]));

        if (distance < tol_sqr)
        {

          auto id_j = udn_2411[j].record1[0];

          similar_points[i].push_back(id_j);
          similar_list.push_back(id_j);
          to_remove.push_back(j);
        }
      }
    }

    auto &udn_2412 = unv_container[uci].udn_2412;
    for (unsigned i = 0; i < udn_2412.size(); ++i)
    {

      int FE_Id = udn_2412[i].record1[1];
      bool beam_type = (FE_Id == 11 || FE_Id == 21 || FE_Id == 22 || FE_Id == 23 || FE_Id == 24);

      unsigned num_of_elements = unv_container[uci].udn_2412[i].record1[5];
      for (unsigned j = 0; j < num_of_elements; ++j)
      {

        for (unsigned m = 0; m < similar_points.size(); ++m)
        {
          for (unsigned n = 0; n < similar_points[m].size(); ++n)
          {

            if (beam_type)
            {
              if (udn_2412[i].record3[j] == similar_points[m][n])
              {
                udn_2412[i].record3[j] = udn_2411[m].record1[0];
              }
            }
            else
            {
              if (udn_2412[i].record2[j] == similar_points[m][n])
              {
                udn_2412[i].record2[j] = udn_2411[m].record1[0];
              }
            }
          }
        }
      }
    }

    std::cout << "no_similar_points: " << to_remove.size() << "\n";

    if (to_remove.size() == 0)
      return;

    std::sort(to_remove.begin(), to_remove.end(), std::greater<unsigned>());
    for (unsigned i = 0; i < to_remove.size(); ++i)
    {
      udn_2411.erase(udn_2411.begin() + to_remove[i]);
    }

    // initial vertex records
    //
    // inx  record1[0]
    // 0:    1
    // 1:    2
    // 2:    3
    // 3:    4
    // ...:  ...
    // n:    n+1

    // AFTER DELETION OF SOME SIMILAR POINTS AND FIXING IDS:

    // inx  old.record1[0]   new.record1[0]
    // 0:    1                1
    // 1:    2                2
    // 2:    5                3
    // 3:    11               4
    // ...:  ...             ...
    // n-m:  o               n-m+1

    // fixing record1[0] from 1 to end
    for (unsigned i = 0; i < udn_2411.size(); ++i)
    {
      int new_id = i + 1;
      new_record1_0[udn_2411[i].record1[0]] = new_id;
      udn_2411[i].record1[0] = new_id;
    }

    for (unsigned i = 0; i < udn_2412.size(); ++i)
    {

      int FE_Id = udn_2412[i].record1[1];
      bool beam_type = (FE_Id == 11 || FE_Id == 21 || FE_Id == 22 || FE_Id == 23 || FE_Id == 24);

      unsigned num_of_elements = unv_container[uci].udn_2412[i].record1[5];
      for (unsigned j = 0; j < num_of_elements; ++j)
      {

        if (beam_type)
        {
          udn_2412[i].record3[j] = new_record1_0[udn_2412[i].record3[j]];
        }
        else
        {
          udn_2412[i].record2[j] = new_record1_0[udn_2412[i].record2[j]];
        }
      }
    }
  }
}
