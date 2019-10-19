
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
#include <array>
#include <string>
#include "vector.h"

int main (int argc, char **argv) {

  const double PI = 3.14159265;

  int start_xyz_frame = 5;

  int atom_type_center = 0;
  int atom_type_target = 1;
  

  double radius_max = 8.0;
  double d_radius = 0.01;
  double d_radius2 = 0.005;
  
  double fix_constant = d_radius2 / d_radius;

  bool has_velocity = false; // the xyz file contain velocity
  
  // PROCESSING ARGUMENTS

  for (int i = 0; i < argc; ++i) {
    //std::cout <<"arg: "<< i << " = " << argv[i] << " v " << (argv[i]=="-v") <<std::endl;
    std::string arg = argv[i];
    if (arg=="-v") {
      has_velocity = true;
      std::cout << "heet" << std::endl;
    }

    if (arg=="-tc") {
      atom_type_center = std::atoi(argv[i+1]);
      ++i;
    }

    if (arg=="-tt") {
      atom_type_target = std::atoi(argv[i+1]);
      ++i;
    }

  }

  if (argc==1)
    std::cout << "no argument from command line. Using default parameters.\n";


  double d = 8.0;
  caviar::Vector<double> small_box_min {-d,-d,-d};
  caviar::Vector<double> small_box_max {+d,+d,+d};

  std::ifstream ifs ("o_xyz.xyz");

  std::ofstream ofs ("o_radial.txt");


  int num_of_radius = radius_max / d_radius;

  std::vector<double> radius_list(num_of_radius, 0.0);
  std::vector<double> density_list(num_of_radius, 0.0);  

  for (int i = 0; i < num_of_radius; ++i) {
    radius_list[i] = d_radius * i;
  }


  // WRITING CODE CONFIGURATIONS
  std::cout << "-- CONFIGURATIONS -- "     << has_velocity << "\n";
  std::cout << "has_velocity: "     << has_velocity << "\n";
  std::cout << "atom_type_center: " << atom_type_center << "\n";
  std::cout << "atom_type_target: " << atom_type_target << "\n";
  std::cout << "\n";
 

  // SYSTEM INITIALIZATONS
  std::vector<int>    type;
  std::vector<caviar::Vector<double> > position;

  std::vector<int> selected_id_list;

  int total_selected_atoms = 0;

  int frame_counter = 0;  

  std::cout << "\n-- PROCESSING -- "     << has_velocity << "\n";
  // PROCESSING THE INPUT FILE
  while (!ifs.eof()) {

    int selected_id = -1;
    selected_id_list.clear();

    int num_of_atoms = 0;
    ifs >> num_of_atoms;

    type.resize(num_of_atoms);
    position.resize(num_of_atoms);



    // don't repeat the last line
    if (ifs.eof()) break; 



    // Ignore the second line. It is contaning "Atom".
    {
      std::string dummyLine;
      getline(ifs, dummyLine);
      getline(ifs, dummyLine);
    } 

   
    // input the frame from file to memory  
    // and select one of the atoms that exist inside our small box. 
    for (int j = 0; j < num_of_atoms; ++j) {
    
      ifs >> type[j];

      double x=0,y=0,z=0, vx=0, vy=0, vz=0;
      if (has_velocity)
        ifs >> x >> y >> z >> vx >> vy >> vz;
      else
        ifs >> x >> y >> z;

      position[j] = caviar::Vector<double> {x, y, z};

      // select the atoms here inside the small box

        if ( small_box_min.x < position[j].x  &&  position[j].x < small_box_max.x
          && small_box_min.y < position[j].y  &&  position[j].y < small_box_max.y 
          && small_box_min.z < position[j].z  &&  position[j].z < small_box_max.z  ) {
        
          if (type[j] == atom_type_center) {
            selected_id_list.push_back(j);
          }
        }
      
   
    }
  
    ++frame_counter;
    if (frame_counter < start_xyz_frame) continue;
    std::cout << "frame number: " << frame_counter << "\n";
    std::cout << "num of selected atoms at this frame:" << selected_id_list.size() << "\n";


    total_selected_atoms += selected_id_list.size();
    // calculating radial distribution from the memory for the atom
    for (int i = 0; i < radius_list.size(); ++i) {

      auto d_v = 4.0*PI*radius_list[i]*radius_list[i]*(2.0*d_radius2);  

      for (int k = 0; k < selected_id_list.size(); ++k) {
        auto center = position [selected_id_list[k]];
        //std::cout << center << std::endl;

        for (int j = 0; j < num_of_atoms; ++j) {

          // XXX: not sure. excluding the selected atom
          if (j == selected_id_list[k]) continue; 

          if (type[j] == atom_type_target) {

            auto dist = position[j] - center;
            auto dist_norm = std::sqrt(dist * dist);

            if (radius_list[i] - d_radius2 <= dist_norm && dist_norm < radius_list[i] + d_radius2) {
              density_list[i] += 1.0 / d_v;
            }

          }
        }
      }
    }



  } 

  std::cout << "\n-- FINISHED --" << "\n";
  std::cout << "total frames:"       << frame_counter << "\n";
  std::cout << "total_selected_atoms:" << total_selected_atoms << "\n";

  for (int i = 0; i < radius_list.size(); ++i) {
    ofs << radius_list[i] << " " << density_list[i]/(total_selected_atoms) << "\n";
  }


}

