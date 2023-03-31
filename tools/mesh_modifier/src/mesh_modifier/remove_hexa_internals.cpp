
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

namespace mesh_modifier
{

  void Mesh_modifier::remove_hexa_internals(const double tol_sqr)
  {
    unsigned uci = unv_container.size() - 1;
    if (uci < 0)
    {
      std::cout << "error: there's no unv_container\n";
      return;
    }
    std::cout << "remove_similar_points: " << std::endl;
    remove_similar_points(tol_sqr);

    std::vector<std::vector<unsigned>> hexa_edges;
    std::vector<std::vector<unsigned>> hexa_quads;

    std::cout << "find_hexahedrons: " << std::endl;
    find_hexahedrons(hexa_edges, hexa_quads);

    auto edges1 = count_edges();
    auto quads1 = count_quads();

    std::vector<unsigned> removed_entity_label;
    std::vector<unsigned> entity_index_to_remove;
    unsigned sum_isolated_entities = 0;
    unsigned sum_internal_entities = 0;
    std::cout << "find isolated entities: " << std::endl;
    auto &udn_2412 = unv_container[uci].udn_2412;
    for (unsigned i = 0; i < udn_2412.size(); ++i)
    {
      int no_entities_found = -1;
      if (udn_2412[i].record1[1] == 11)
      {
        std::vector<unsigned> e{udn_2412[i].record3[0], udn_2412[i].record3[1]};
        if (e[0] > e[1])
        {
          unsigned tmp = e[0];
          e[0] = e[1];
          e[1] = tmp;
        }

        no_entities_found = std::count(hexa_edges.begin(), hexa_edges.end(), e);
      }
      else if (udn_2412[i].record1[1] == 44)
      {
        std::vector<unsigned> e{udn_2412[i].record2[0], udn_2412[i].record2[1], udn_2412[i].record2[2], udn_2412[i].record2[3]};

        std::sort(e.begin(), e.end());
        no_entities_found = std::count(hexa_quads.begin(), hexa_quads.end(), e);
      }

      if (no_entities_found == 0)
      {
        entity_index_to_remove.push_back(i);
        removed_entity_label.push_back(udn_2412[i].record1[0]);
        ++sum_isolated_entities;
      }
      if (no_entities_found > 1)
      {
        entity_index_to_remove.push_back(i);
        removed_entity_label.push_back(udn_2412[i].record1[0]);
        ++sum_internal_entities;
      }
    }

    std::sort(entity_index_to_remove.begin(), entity_index_to_remove.end(), std::greater<unsigned>());
    for (unsigned i = 0; i < entity_index_to_remove.size(); ++i)
      udn_2412.erase(udn_2412.begin() + entity_index_to_remove[i]);

    std::cout << "Number of removed isolated entities: "
              << sum_isolated_entities << "\n";
    std::cout << "Number of removed internal entities: "
              << sum_internal_entities << "\n";

    /*

    std::cout << "remove_internal_edges: " << std::endl;
    std::vector<unsigned> removed_entity_label;
    for (unsigned i = 0 ; i < hexa_edges.size(); ++i) {
      unsigned no_edge =  std::count (hexa_edges.begin(),hexa_edges.end(), hexa_edges[i]);
      if (no_edge > 1) {
        remove_edge (hexa_edges[i], removed_entity_label);
      }
    }


    std::cout << "hexa_quads.size(): " << hexa_quads.size() << std::endl;
    std::cout << "remove_quad: " << std::endl;
    for (unsigned i = 0 ; i < hexa_quads.size(); ++i) {
      unsigned no_quad =  std::count (hexa_quads.begin(),hexa_quads.end(), hexa_quads[i]);
      //std::cout << " i : " << i << " no_quads: " << no_quad << "\n";
      if (no_quad > 1) {
        remove_quad (hexa_quads[i], removed_entity_label);
      }
    }*/

    auto edges2 = count_edges();
    auto quads2 = count_quads();

    std::cout << "number of deleted edges: " << edges1 - edges2 << "\n";
    std::cout << "number of deleted quads: " << quads1 - quads2 << "\n";

    std::cout << "remove_hexa_internal_groups: " << std::endl;
    remove_hexa_internal_groups(removed_entity_label);
  }

  void Mesh_modifier::remove_hexa_internal_groups(std::vector<unsigned> &removed_entity_label)
  {
    unsigned uci = unv_container.size() - 1;
    auto &udn_2467 = unv_container[uci].udn_2467;
    auto &removed_labels = removed_entity_label;
    std::cout << "removed_labels: " << removed_labels.size() << "\n";
    std::cout << "udn_2467.size(): " << udn_2467.size() << "\n";
    for (int i = udn_2467.size() - 1; i > -1; --i)
    {
      std::cout << i << "\n";
      unsigned num_of_entities = udn_2467[i].record1[7];
      for (int j = num_of_entities - 1; j > -1; --j)
      {
        unsigned m = (4 * j) + 1;
        unsigned tag = udn_2467[i].record3[m];

        for (unsigned k = 0; k < removed_labels.size(); ++k)
        {
          if (tag == removed_labels[k])
          {
            udn_2467[i].record3.erase(udn_2467[i].record3.begin() + (4 * j) + 3);
            udn_2467[i].record3.erase(udn_2467[i].record3.begin() + (4 * j) + 2);
            udn_2467[i].record3.erase(udn_2467[i].record3.begin() + (4 * j) + 1);
            udn_2467[i].record3.erase(udn_2467[i].record3.begin() + (4 * j) + 0);
            --udn_2467[i].record1[7];
            std::cout << "Warning: group (udn_2467) number " << udn_2467[i].record1[0]
                      << " with name: " << udn_2467[i].record2
                      << " element n. " << j << " has been removed due being internal.\n";
            break;
          }
        }
      }

      if (udn_2467[i].record1[7] == 0)
      {
        std::cout << "Warning: group (udn_2467) number " << udn_2467[i].record1[0]
                  << " with name " << udn_2467[i].record2
                  << " removed completely, due being internal.\n";
        udn_2467.erase(udn_2467.begin() + i);
      }
    }

    removed_labels.clear();
  }

  unsigned Mesh_modifier::count_edges()
  {
    unsigned uci = unv_container.size() - 1;
    auto &udn_2412 = unv_container[uci].udn_2412;
    unsigned sum = 0;
    for (unsigned i = 0; i < udn_2412.size(); ++i)
      if (udn_2412[i].record1[1] == 11)
        ++sum;
    return sum;
  }

  unsigned Mesh_modifier::count_quads()
  {
    unsigned uci = unv_container.size() - 1;
    auto &udn_2412 = unv_container[uci].udn_2412;
    unsigned sum = 0;
    for (unsigned i = 0; i < udn_2412.size(); ++i)
      if (udn_2412[i].record1[1] == 44)
        ++sum;

    return sum;
  }

  void Mesh_modifier::find_hexahedrons(std::vector<std::vector<unsigned>> &hexa_edges, std::vector<std::vector<unsigned>> &hexa_quads)
  {
    unsigned uci = unv_container.size() - 1;
    unsigned num_of_2412 = unv_container[uci].udn_2412.size();
    auto &udn_2411 = unv_container[uci].udn_2411;
    for (unsigned i = 0; i < num_of_2412; ++i)
    {
      auto &udn_2412 = unv_container[uci].udn_2412[i];
      if (udn_2412.record1[1] == 115)
      {
        auto &pl = udn_2412.record2;

        add_edge_to_vector(hexa_edges, pl[0], pl[1]);
        add_edge_to_vector(hexa_edges, pl[1], pl[2]);
        add_edge_to_vector(hexa_edges, pl[2], pl[3]);
        add_edge_to_vector(hexa_edges, pl[3], pl[0]);

        add_edge_to_vector(hexa_edges, pl[4], pl[5]);
        add_edge_to_vector(hexa_edges, pl[5], pl[6]);
        add_edge_to_vector(hexa_edges, pl[6], pl[7]);
        add_edge_to_vector(hexa_edges, pl[7], pl[4]);

        add_edge_to_vector(hexa_edges, pl[0], pl[4]);
        add_edge_to_vector(hexa_edges, pl[1], pl[5]);
        add_edge_to_vector(hexa_edges, pl[2], pl[6]);
        add_edge_to_vector(hexa_edges, pl[3], pl[7]);

        add_quad_to_vector(hexa_quads, pl[0], pl[1], pl[2], pl[3]);
        add_quad_to_vector(hexa_quads, pl[4], pl[5], pl[6], pl[7]);
        add_quad_to_vector(hexa_quads, pl[0], pl[1], pl[5], pl[4]);
        add_quad_to_vector(hexa_quads, pl[2], pl[6], pl[7], pl[3]);
        add_quad_to_vector(hexa_quads, pl[1], pl[5], pl[6], pl[2]);
        add_quad_to_vector(hexa_quads, pl[0], pl[4], pl[7], pl[3]);
      }
    }
  }

}
