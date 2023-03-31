
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
#include <algorithm>
#include <iterator>
#include <iomanip>

namespace mesh_modifier
{

  void Mesh_modifier::import_vtk(const std::string &filename, bool)
  {

    unv_container.push_back(Unv_container());
    unsigned uci = unv_container.size() - 1;

    std::ifstream ifs;
    ifs.open(filename.c_str());

    std::string c;

    while (true)
    {
      ifs >> c;
      std::cout << c << "\n";
      if (ifs.eof())
        break;
      if (c == "POINTS")
        break;
    }
    /*
        ifs >> c;

        if (c == "#")
          ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if (c == "ASCII")
          ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if (c == "DATASET") // TODO: add the other cases
          ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        */
    if (c == "POINTS")
    {
      unsigned int num_points;
      ifs >> num_points;
      ifs >> c; // the type of points such as 'float'. Does not matter.

      for (unsigned int i = 0; i < num_points; ++i)
      {
        double x, y, z;

        ifs >> x;
        ifs >> y;
        ifs >> z;

        unv_container[uci].udn_2411.push_back(Universal_dataset_number_2411());
        unsigned udi = unv_container[uci].udn_2411.size() - 1;

        unv_container[uci].udn_2411[udi].record1[0] = i + 1;
        unv_container[uci].udn_2411[udi].record1[1] = 1;
        unv_container[uci].udn_2411[udi].record1[2] = 1;
        unv_container[uci].udn_2411[udi].record1[3] = 11;

        unv_container[uci].udn_2411[udi].record2[0] = x;
        unv_container[uci].udn_2411[udi].record2[1] = y;
        unv_container[uci].udn_2411[udi].record2[2] = z;
        std::cout << x << " " << y << " " << z << "\n";
      }
    }

    ifs >> c;

    if (c == "POLYGONS")
    {

      unsigned int num_of_poly;
      ifs >> num_of_poly;

      unsigned int num_of_num;
      ifs >> num_of_num;
      std::cout << "\nPOLYGONS " << num_of_poly << " " << num_of_num << "\n";
      for (unsigned int i = 0; i < num_of_poly; ++i)
      {

        unv_container[uci].udn_2412.push_back(Universal_dataset_number_2412());
        unsigned udi = unv_container[uci].udn_2412.size() - 1;

        unv_container[uci].udn_2412[udi].record1.push_back(i + 1);
        unv_container[uci].udn_2412[udi].record1.push_back(41); // only 41 is supported
        unv_container[uci].udn_2412[udi].record1.push_back(2);
        unv_container[uci].udn_2412[udi].record1.push_back(1);
        unv_container[uci].udn_2412[udi].record1.push_back(7);
        unv_container[uci].udn_2412[udi].record1.push_back(3);
        unsigned int tmp;
        unsigned int p1, p2, p3;
        ifs >> tmp; // TODO: complete if 'tmp != 3'

        ifs >> p1;
        ifs >> p2;
        ifs >> p3;

        unv_container[uci].udn_2412[udi].record2.push_back(p1 + 1);
        unv_container[uci].udn_2412[udi].record2.push_back(p2 + 1);
        unv_container[uci].udn_2412[udi].record2.push_back(p3 + 1);
        std::cout << tmp << " " << p1 << " " << p2 << " " << p3 << "\n";
      }
    }

    ifs.close();
  }

}
