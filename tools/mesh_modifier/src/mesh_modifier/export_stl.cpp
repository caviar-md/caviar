
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
#include <cmath>

namespace mesh_modifier
{

  void Mesh_modifier::export_stl_boundary(const std::string &filename)
  {

    unsigned uci = unv_container.size() - 1;
    if (uci < 0)
    {
      std::cout << "error: there's no unv_container\n";
      return;
    }

    std::ofstream ofs;
    ofs.open(filename.c_str());

    export_stl_headers(ofs);
    export_stl_polygons(ofs);
    export_stl_enders(ofs);

    ofs.close();
  }

  void Mesh_modifier::export_stl_headers(std::ofstream &ofs)
  {
    ofs << "solid shape, STL ascii file, created with Mesh_modifier code\n";
  }

  void Mesh_modifier::export_stl_enders(std::ofstream &ofs)
  {
    ofs << "endsolid shape";
  }

  void Mesh_modifier::export_stl_polygons(std::ofstream &ofs)
  {
    unsigned uci = unv_container.size() - 1;
    auto &udn_2411 = unv_container[uci].udn_2411;
    auto &udn_2412 = unv_container[uci].udn_2412;

    ofs << std::scientific << std::setprecision(6);

    for (int i = 0; i < udn_2412.size(); ++i)
    {

      // plane triangle cases
      if (udn_2412[i].record1[1] == 41)
      {
        auto i1 = unv_container[uci].udn_2412[i].record2[0];
        auto i2 = unv_container[uci].udn_2412[i].record2[1];
        auto i3 = unv_container[uci].udn_2412[i].record2[2];

        auto v1 = udn_2411[i1].record2;
        auto v2 = udn_2411[i2].record2;
        auto v3 = udn_2411[i3].record2;
        Vector<double> vv1{v1[0], v1[1], v1[2]};
        Vector<double> vv2{v2[0], v2[1], v2[2]};
        Vector<double> vv3{v3[0], v3[1], v3[2]};
        auto normal = cross_product(vv2 - vv1, vv3 - vv1);
        normal = normal * (1.0 / std::sqrt(normal * normal));

        ofs << " facet normal " << normal.x << " " << normal.y << " " << normal.z << "\n";
        ofs << "   outer loop\n";
        ofs << "     vertex " << vv1.x << " " << vv1.y << " " << vv1.z << "\n";
        ofs << "     vertex " << vv2.x << " " << vv2.y << " " << vv2.z << "\n";
        ofs << "     vertex " << vv3.x << " " << vv3.y << " " << vv3.z << "\n";
        ofs << "   endloop\n";
        ofs << " endfacet\n";

        // plane Linear Quadrilateral
      }
      else if (udn_2412[i].record1[1] == 44)
      {
        auto i1 = unv_container[uci].udn_2412[i].record2[0] - 1;
        auto i2 = unv_container[uci].udn_2412[i].record2[1] - 1;
        auto i3 = unv_container[uci].udn_2412[i].record2[2] - 1;
        auto i4 = unv_container[uci].udn_2412[i].record2[3] - 1;

        auto v1 = udn_2411[i1].record2;
        auto v2 = udn_2411[i2].record2;
        auto v3 = udn_2411[i3].record2;
        auto v4 = udn_2411[i4].record2;
        Vector<double> vv1{v1[0], v1[1], v1[2]};
        Vector<double> vv2{v2[0], v2[1], v2[2]};
        Vector<double> vv3{v3[0], v3[1], v3[2]};
        Vector<double> vv4{v4[0], v4[1], v4[2]};

        auto normal = cross_product(vv2 - vv1, vv3 - vv1);
        normal = normal * (1.0 / std::sqrt(normal * normal));

        ofs << " facet normal " << normal.x << " " << normal.y << " " << normal.z << "\n";
        ofs << "   outer loop\n";
        ofs << "     vertex " << vv1.x << " " << vv1.y << " " << vv1.z << "\n";
        ofs << "     vertex " << vv2.x << " " << vv2.y << " " << vv2.z << "\n";
        ofs << "     vertex " << vv3.x << " " << vv3.y << " " << vv3.z << "\n";
        ofs << "   endloop\n";
        ofs << " endfacet\n";

        ofs << " facet normal " << normal.x << " " << normal.y << " " << normal.z << "\n";
        ofs << "   outer loop\n";
        ofs << "     vertex " << vv1.x << " " << vv1.y << " " << vv1.z << "\n";
        ofs << "     vertex " << vv3.x << " " << vv3.y << " " << vv3.z << "\n";
        ofs << "     vertex " << vv4.x << " " << vv4.y << " " << vv4.z << "\n";
        ofs << "   endloop\n";
        ofs << " endfacet\n";
      }
      else
      {
        // not supported yet!
      }
    }
  }

} // mesh_modifier
