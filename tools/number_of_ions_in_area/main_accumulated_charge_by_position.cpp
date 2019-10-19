
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

int main (int argc, char **argv) {
  bool has_velocity = false; // the xyz file contain velocity
  int xyz_steps = 10000;
  double time_converter = 0.000006; // step to time in (ps)


  //std::cout << "argc: " << argc << "\nargv:";
  for (int i = 0; i < argc; ++i) {
    std::cout << " " << argv[i] << "\n";

    std::string arg = argv[i];
    if (arg=="-v") {
      has_velocity = true;
    }

    if (arg=="-s") {
      if (i+1 >= argc) {
        std::cout << " Error: expected int value for 'xyz_steps'\n";
        return 1;
      }
      xyz_steps = std::atoi(argv[i+1]);
      ++i;
    }
    if (arg=="-t") {
      if (i+1 >= argc) {
        std::cout << " Error: expected double value for 'time_converter'\n";
        return 1;
      }
      time_converter = std::atof(argv[i+1]);
      ++i;
    }
  }

  if (argc==1)
    std::cout << "no argument from input file. Using default parameters.\n";

  std::cout << "has_velocity: " << has_velocity << "\n";
  std::cout << "xyz_steps: " << xyz_steps << "\n";
  std::cout << "time_converter: " << time_converter << "\n";

  bool output_total_number_of_ions = true;
  int max_atom_type = 2;
  int start_process = 0;
  //double dt = 0.001; // in simulation units



  std::ifstream ifs ("o_xyz.xyz");
  std::ofstream ofs ("o_accumulated_charge");

  //std::vector<double> sum_values(num_of_columns - 1, 0);

  int frame_counter = 0;  
  std::vector < std::array <double,2> > x_positions;
  x_positions.push_back({-31,0});
  x_positions.push_back({-31,-16});
  x_positions.push_back({-16,-8});
  x_positions.push_back({-8,0});
  x_positions.push_back({0,24});


  std::cout << "calculating number of ions at the x intervals: ";
  for (int i = 0; i < x_positions.size(); ++i) {
    std::cout << " [ " << x_positions[i][0] <<" : " <<x_positions[i][1] << " ] " << ","; 
  }

  std::cout << std::endl;
  std::vector<std::vector <int> > x_number;
  x_number.resize(x_positions.size());
  for (auto && i : x_number)
    i.resize(max_atom_type);

  std::ofstream ofs_t0("o_t0_ionnumber");  
  std::ofstream ofs_t1("o_t1_ionnumber");

  std::ofstream ofs_sum;
  if (output_total_number_of_ions)
    ofs_sum.open("o_sum_ionnumber");

  while (!ifs.eof()) {

    for (int k = 0; k < x_positions.size()-1; ++k) {
      x_number[k][0] = 0;
      x_number[k][1] = 0;
    }


    int num_of_atoms = 0;
    ifs >> num_of_atoms;
   
    if (ifs.eof()) break; // don't repeat the last line

    {
      // Ignore the second line    
      std::string dummyLine;
      getline(ifs, dummyLine);
      getline(ifs, dummyLine);
    } 

   
   
    for (int j = 0; j < num_of_atoms; ++j) {
    
      int type;
      ifs >> type;
      if (type+1 > max_atom_type) {
        std::cout << "max atom type is larger than " << max_atom_type << std::endl;
        //return 1;
      }
       //std::cout << type << " ";

      double x=0,y=0,z=0, vx=0, vy=0, vz=0;
      if (has_velocity)
        ifs >> x >> y >> z >> vx >> vy >> vz;
      else
        ifs >> x >> y >> z;

       //std::cout << x << " " << y << " " << z <<  "\n";      

      for (int k = 0; k < x_positions.size(); ++k) {
        if (x_positions[k][0] < x && x < x_positions[k][1]) {
          x_number[k][type] ++;
        }                
      }

    
    }  

    double time = time_converter * frame_counter * xyz_steps;


    ofs_t0 << time ;
    ofs_t1 << time ;

    int sum_t0 = 0;
    int sum_t1 = 0;

    for (int k = 0; k < x_number.size(); ++k) {
      sum_t0 +=x_number[k][0];
      sum_t1 +=x_number[k][1];
    }

    if (!(sum_t0==0 && sum_t1==0)) {

      for (int k = 0; k < x_number.size(); ++k) {
        ofs_t0 << " " << x_number[k][0];
        ofs_t1 << " " << x_number[k][1];

      }
      ofs_t0 << " " << sum_t0 << "\n" ;
      ofs_t1 << " " << sum_t1 << "\n" ;
    }


    ++frame_counter;
  } 



}

