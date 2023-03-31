
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

#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <iomanip>

namespace mesh_modifier
{

  void Mesh_modifier::export_vtk_boundary(const std::string &filename)
  {

    unsigned uci = unv_container.size() - 1;
    if (uci < 0)
    {
      std::cout << "error: there's no unv_container\n";
      return;
    }

    std::ofstream ofs;
    ofs.open(filename.c_str());

    export_vtk_headers(ofs);
    export_vtk_points(ofs);
    export_vtk_polygons(ofs);

    ofs.close();
  }

  void Mesh_modifier::export_vtk_polygons(std::ofstream &ofs)
  {
    unsigned uci = unv_container.size() - 1;
    auto &udn_2411 = unv_container[uci].udn_2411;
    auto &udn_2412 = unv_container[uci].udn_2412;

    int no_41 = 0, no_44 = 0;

    for (int i = 0; i < udn_2412.size(); ++i)
    {
      if (udn_2412[i].record1[1] == 41)
        ++no_41;
      else if (udn_2412[i].record1[1] == 44)
        ++no_44;
    }
    int tot_polygons = no_41 + 2 * no_44;

    ofs << "\nCELLS " << tot_polygons << " " << 4 * tot_polygons << "\n";

    for (int i = 0; i < udn_2412.size(); ++i)
    {

      // plane triangle cases
      if (udn_2412[i].record1[1] == 41)
      {
        auto i1 = unv_container[uci].udn_2412[i].record2[0] - 1;
        auto i2 = unv_container[uci].udn_2412[i].record2[1] - 1;
        auto i3 = unv_container[uci].udn_2412[i].record2[2] - 1;

        ofs << "3 " << i1 << " " << i2 << " " << i3 << "\n";
        // plane Linear Quadrilateral
      }
      else if (udn_2412[i].record1[1] == 44)
      {
        auto i1 = unv_container[uci].udn_2412[i].record2[0] - 1;
        auto i2 = unv_container[uci].udn_2412[i].record2[1] - 1;
        auto i3 = unv_container[uci].udn_2412[i].record2[2] - 1;
        auto i4 = unv_container[uci].udn_2412[i].record2[3] - 1;

        ofs << "3 " << i1 << " " << i2 << " " << i3 << "\n";
        ofs << "3 " << i1 << " " << i3 << " " << i4 << "\n";
      }
      else
      {
        // not supported yet!
      }
    }
  }

  void Mesh_modifier::export_vtk(const std::string &filename)
  {

    unsigned uci = unv_container.size() - 1;
    if (uci < 0)
    {
      std::cout << "error: there's no unv_container\n";
      return;
    }

    std::ofstream ofs;
    ofs.open(filename.c_str());

    export_vtk_headers(ofs);
    export_vtk_points(ofs);
    export_vtk_cells(ofs);
    export_vtk_cell_types(ofs);
    export_vtk_cell_data(ofs);

    ofs.close();
  }

  void Mesh_modifier::export_vtk_headers(std::ofstream &ofs)
  {
    std::string text[4];
    text[0] = "# vtk DataFile Version 4.1\n";
    text[1] = "vtk output\n";
    text[2] = "ASCII\n";
    text[3] = "DATASET POLYDATA\n";
    // text[3] = "DATASET UNSTRUCTURED_GRID\n";

    for (unsigned int i = 0; i < 4; ++i)
    {
      ofs << text[i];
    }
  }

  void Mesh_modifier::export_vtk_points(std::ofstream &ofs)
  {
    unsigned uci = unv_container.size() - 1;
    unsigned no_points = unv_container[uci].udn_2411.size();
    // ofs << "POINTS " << no_points << " double\n";
    ofs << "POINTS " << no_points << " float\n";
    for (unsigned i = 0; i < unv_container[uci].udn_2411.size(); ++i)
    {
      ofs << unv_container[uci].udn_2411[i].record2[0] << " "
          << unv_container[uci].udn_2411[i].record2[1] << " "
          << unv_container[uci].udn_2411[i].record2[2] << "\n";
    }
  }

  // http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html
  // The CELLS keyword requires two parameters: the number of cells n and the
  // size of the cell list size. The cell list size is the total number of
  // integer values required to represent the list (i. e., sum
  // of numPoints and connectivity indices over each cell)
  /*
  http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html
  DATASET UNSTRUCTURED_ GRID
  POINTS n dataType
  p 0x p 0y p 0z
  p 1x p 1y p 1z

  ...
  p (n-1) x p (n-1) y p (n-1) z

  CELLS n size
  numPoints 0 ,i, j, k, l,...
  numPoints 1 ,i, j, k, l,...
  numPoints 2 ,i, j, k, l,...
  ...
  numPoints n-1 ,i, j, k, l,...
  Simple Legacy Formats 5
  CELL_ TYPES n
  type 0
  type 1
  type 2
  */
  void Mesh_modifier::export_vtk_cells(std::ofstream &ofs)
  {
    unsigned uci = unv_container.size() - 1;
    auto &udn_2412 = unv_container[uci].udn_2412;
    unsigned no_cells = 0, no_ints = 0;

    for (auto i : udn_2412)
    {
      auto cell_type = i.record1[1];
      if (cell_type == 41)
      {
        ++no_cells;
        no_ints += 4;
      }
      if (cell_type == 44)
      {
        ++no_cells;
        no_ints += 5;
      }
      if (cell_type == 115)
      {
        ++no_cells;
        no_ints += 9;
      }
    }

    // ofs << "\nCELLS " << no_cells << " " << no_ints << "\n";
    ofs << "\nPOLYGONS " << no_cells << " " << no_ints << "\n";

    std::vector<unsigned> li;
    make_label_to_index(unv_container[uci].udn_2411, li);

    for (unsigned i = 0; i < udn_2412.size(); ++i)
    {
      auto cell_type = udn_2412[i].record1[1];
      if (cell_type == 115)
      {
        ofs << "8";
        for (unsigned int j = 0; j < 8; ++j)
          ofs << " " << li[udn_2412[i].record2[j]];
        ofs << "\n";
      }
    }

    for (unsigned i = 0; i < udn_2412.size(); ++i)
    {
      auto cell_type = udn_2412[i].record1[1];
      if (cell_type == 41)
      {
        ofs << "3";
        for (unsigned int j = 0; j < 3; ++j)
          ofs << " " << li[udn_2412[i].record2[j]];
        ofs << "\n";
      }
    }

    for (unsigned i = 0; i < udn_2412.size(); ++i)
    {
      auto cell_type = udn_2412[i].record1[1];
      if (cell_type == 44)
      {
        ofs << "4";
        for (unsigned int j = 0; j < 4; ++j)
          ofs << " " << li[udn_2412[i].record2[j]];
        ofs << "\n";
      }
    }
  }

  void Mesh_modifier::export_vtk_cell_types(std::ofstream &ofs)
  {
    unsigned uci = unv_container.size() - 1;
    auto &udn_2412 = unv_container[uci].udn_2412;
    unsigned no_cells = 0;

    for (auto i : udn_2412)
    {
      auto cell_type = i.record1[1];
      if (cell_type == 44 || cell_type == 115)
      {
        ++no_cells;
      }
    }

    ofs << "\nCELL_TYPES " << no_cells << "\n ";

    for (auto i : udn_2412)
    {
      auto cell_type = i.record1[1];
      if (cell_type == 115)
      {
        ofs << "12 ";
      }
    }
    for (auto i : udn_2412)
    {
      auto cell_type = i.record1[1];
      if (cell_type == 44)
      {
        ofs << "9 ";
      }
    }
  }

  void Mesh_modifier::export_vtk_cell_data(std::ofstream &ofs)
  {
    unsigned uci = unv_container.size() - 1;
    auto &udn_2412 = unv_container[uci].udn_2412;
    auto &udn_2467 = unv_container[uci].udn_2467;
    unsigned no_cells = udn_2412.size();
    if (udn_2467.size() == 0)
      return;

    std::vector<unsigned> li_cell, li_quad;
    unsigned max_label = 0;
    for (unsigned int i = 0; i < udn_2412.size(); ++i)
    {
      if ((udn_2412[i].record1[1]) == 44)
        li_quad.push_back(i);
      else if ((udn_2412[i].record1[1]) == 115)
        li_cell.push_back(i);
    }

    auto no_data = li_quad.size() + li_cell.size();
    ofs << "\nCELL_DATA " << no_data << "\n";
    ofs << "SCALARS MaterialID double\n";
    ofs << "LOOKUP_TABLE default\n";
    /*
        std::vector<unsigned> li;
        make_label_to_index (unv_container[uci].udn_2412, li);

        std::vector<unsigned> ids (udn_2412.size(), 0);
        for (unsigned int i = 0; i < udn_2467.size(); ++i) {
          unsigned num_of_elements = udn_2467[i].record1[7];
          for (unsigned j = 0; j < num_of_elements; ++j) {
            unsigned m = 4*j + 1;
            unsigned label = udn_2467[i].record3[m];
            ids [li[label]] = udn_2467[j].record2;
          }
        }
    */
    std::vector<unsigned> ids(udn_2412.size(), 0);
    extract_udn_2467_ids(ids);

    for (auto i : li_cell)
      ofs << ids[i] << " ";
    for (auto i : li_quad)
      ofs << ids[i] << " ";
    ofs << "\n";
  }

  void Mesh_modifier::extract_udn_2467_ids(std::vector<unsigned> &ids)
  {

    unsigned uci = unv_container.size() - 1;
    unsigned num_of_2467 = unv_container[uci].udn_2467.size();
    auto &udn_2467 = unv_container[uci].udn_2467;

    std::vector<unsigned> li;
    make_label_to_index(unv_container[uci].udn_2412, li);

    for (unsigned i = 0; i < num_of_2467; ++i)
    {
      unsigned num_of_elements = udn_2467[i].record1[7];
      for (unsigned j = 0; j < num_of_elements; ++j)
      {
        unsigned m = 4 * j + 1;
        unsigned label = udn_2467[i].record3[m];
        ids[li[label]] = udn_2467[j].record2;
      }
    }
  }

}
