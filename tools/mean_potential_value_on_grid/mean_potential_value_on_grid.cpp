
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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// this program calculates mean of induced charges
// read the comments and you will know what to do.
int main()
{

  std::string input_file = "output_pot_3";

  std::string output_profile_file = "potential_profile";

  // std::string output_data_file = "potential_data";

  char profile_direction = 'x'; // x, y or z

  int potential_column = 12;

  //------------
  int frame_col = 1;

  int point_index_col = 2;

  int x_index_col = 3;

  int y_index_col = 4;

  int z_index_col = 5;

  int x_col = 6;

  int y_col = 7;

  int z_col = 8;

  //--------------------------------------------------

  std::ifstream ifs;

  ifs.open(input_file.c_str());

  int default_grid_index = 1000;

  bool set_first_frame_number = true;
  bool in_first_frame = true;
  int first_frame_number = -1;

  std::vector<double> x_vector(default_grid_index, 0.0),
      y_vector(default_grid_index, 0.0),
      z_vector(default_grid_index, 0.0);

  std::vector<double> sum_potential(default_grid_index, 0.0);

  int max_grid_index_x = 1;
  int max_grid_index_y = 1;
  int max_grid_index_z = 1;

  // ---------------------------------------------------------------------
  // importing loop
  // ---------------------------------------------------------------------
  int frame_counter = 0;
  int old_frame = -1;

  int num_of_points_in_one_cross_section = 0;
  int first_index_in_cross_section = -1;
  bool set_first_index_in_cross_section = true;

  while (!ifs.eof())
  // for (int k = 0; k < 200; ++k)
  {

    int frame, point_index, x_index, y_index, z_index;
    double x, y, z;

    ifs >> frame;

    std::cout << frame << " ";
    if (ifs.eof())
      break; // break if eof is reached

    if (old_frame != frame)
    {
      old_frame = frame;
      frame_counter++;
    }

    ifs >> point_index >> x_index >> y_index >> z_index >> x >> y >> z;

    // std::cout << x_index  << " " << y_index << " " <<z_index << " " << x <<" " << y <<" " << z <<" ";

    if (set_first_frame_number)
    {
      first_frame_number = frame;
      set_first_frame_number = false;
    }

    if (frame != first_frame_number)
      in_first_frame = false;

    if (in_first_frame)
    {
      if (max_grid_index_x < x_index)
        max_grid_index_x = x_index;
      if (max_grid_index_y < y_index)
        max_grid_index_y = y_index;
      if (max_grid_index_z < z_index)
        max_grid_index_z = z_index;

      if (x_vector.size() < x_index)
        x_vector.resize(x_index, 0.0);
      if (y_vector.size() < y_index)
        y_vector.resize(y_index, 0.0);
      if (z_vector.size() < z_index)
        z_vector.resize(z_index, 0.0);

      switch (profile_direction)
      {
      case 'x':
        if (sum_potential.size() < x_index)
          sum_potential.resize(x_index, 0.0);
        break;
      case 'y':
        if (sum_potential.size() < y_index)
          sum_potential.resize(y_index, 0.0);
        break;
      case 'z':
        if (sum_potential.size() < z_index)
          sum_potential.resize(z_index, 0.0);
        break;
      default:
        break;
      }

      x_vector[x_index] = x;
      y_vector[y_index] = y;
      z_vector[z_index] = z;

      ///

      if (set_first_index_in_cross_section)
      {
        switch (profile_direction)
        {
        case 'x':
          first_index_in_cross_section = x_index;
          break;
        case 'y':
          first_index_in_cross_section = y_index;
          break;
        case 'z':
          first_index_in_cross_section = z_index;
          break;
        default:
          break;
        }
        set_first_index_in_cross_section = false;
        num_of_points_in_one_cross_section++;
      }
      else
      {
        switch (profile_direction)
        {
        case 'x':
          if (first_index_in_cross_section == x_index)
            num_of_points_in_one_cross_section++;
          break;
        case 'y':
          if (first_index_in_cross_section == y_index)
            num_of_points_in_one_cross_section++;
          break;
        case 'z':
          if (first_index_in_cross_section == z_index)
            num_of_points_in_one_cross_section++;
          break;
        default:
          break;
        }
      }
    }

    double tmp = 0;
    for (int i = 9; i < potential_column + 1; ++i)
    {
      ifs >> tmp;
      std::cout << tmp << " ";
    }

    switch (profile_direction)
    {
    case 'x':
      sum_potential[x_index] += tmp;
      break;
    case 'y':
      sum_potential[y_index] += tmp;
      break;
    case 'z':
      sum_potential[z_index] += tmp;
      break;
    default:
      break;
    }

    std::cout << tmp << std::endl;

    std::cout << " x_index " << x_index << std::endl;
    std::cout << " max_grid_index_x " << max_grid_index_x << std::endl;
  }

  ifs.close();

  // -------------------------------------------------------------------
  // process and export
  //--------------------------------------------------------------------

  // shrinking the size of 'sum_potential'
  switch (profile_direction)
  {
  case 'x':
    sum_potential.resize(max_grid_index_x, 0.0);
    break;
  case 'y':
    sum_potential.resize(max_grid_index_y, 0.0);
    break;
  case 'z':
    sum_potential.resize(max_grid_index_z, 0.0);
    break;
  default:
    break;
  }

  std::ofstream ofs;

  ofs.open(output_profile_file.c_str());

  std::cout << "num_of_points_in_one_cross_section " << num_of_points_in_one_cross_section << std::endl;
  std::cout << "frame_counter " << frame_counter << std::endl;
  int total_num_of_points_in_one_cross_section = num_of_points_in_one_cross_section * frame_counter;
  for (int i = 0; i < sum_potential.size(); ++i)
  {

    switch (profile_direction)
    {
    case 'x':
      ofs << x_vector[i] << " " << sum_potential[i] / total_num_of_points_in_one_cross_section << "\n";
      break;
    case 'y':
      ofs << y_vector[i] << " " << sum_potential[i] / total_num_of_points_in_one_cross_section << "\n";
      break;
    case 'z':
      ofs << z_vector[i] << " " << sum_potential[i] / total_num_of_points_in_one_cross_section << "\n";
      break;
    default:
      break;
    }
  }

  ofs.close();
}
