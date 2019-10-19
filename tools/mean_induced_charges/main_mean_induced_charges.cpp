
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
int main () {

  // ----------------------------------------------------------
  // User variables
  // ----------------------------------------------------------

  // this is a conversion coefficient optional for user. Since the simulation
  // units is usually different from report or paper units, one can convert
  // it using this variable.
  const double conversion_coef_charge = 1.0; 
  const double conversion_coef_time = 1.0; 


  // if true, the converted values of the inputted file will be generated.
  // There's another file consist of the value of induced charges per unit area.
  bool output_converted = true;


  // One may want to start the calculations after some equilibrium time.
  // by setting one of these two variables, the code will follow your lead
  // on when to start the calculations/
  // if it is equal to '-1', it is not activated
  // if both of them are '-1', the calculations starts from the first timestep
  const int start_at_timestep = -1;   
  const int start_at_line = -1;       


  // Number of columns depends on the problem you are simulating.
  // the first colomn is timestep, the rest are induced charges on the areas
  // by ids. The last two columns are sum of charges and sum of abs of charges.
  // this should be set depending on the number of columns of the input file,
  // i.e. 'o_induced_charge'
  const unsigned num_of_columns = 10;       


  // input file name for induced charge.
  const std::string induced_charge = "o_induced_charge"; 


  // input file name for mesh boundary area.
  const std::string mesh_boundary_area = "o_mesh_boundary_area"; 
  

  // output file names for differential capacitance. 
  const std::string mean_induced_charge = "o_mean_induced_charge"; 

  // output file names for differential capacitance. 
  const std::string mean_induced_charge_per_area = "o_mean_induced_charge_per_area"; 

  // output file names for converted induced charge
  const std::string mean_converted_induced_charge = "o_converted_induced_charge"; 

  // output file names for converted induced charge per area. 
  const std::string mean_converted_induced_charge_per_area = "o_converted_induced_charge_per_area"; 

  // ----------------------------------------------------------
  // Calculations starts from here. 
  // ----------------------------------------------------------






  std::ifstream ifs_area (mesh_boundary_area.c_str()); 

  // importing area of the boundaries by id
  std::vector<int> area_id (num_of_columns+2,0);

  // store the inverse of the area
  std::vector<double> area_inv (num_of_columns+2,0);
  double sum_area = 0.0;
  for (unsigned int i =0 ; i < num_of_columns; ++i) {
    int id;
    double a;
    ifs_area >> id;    
    ifs_area >> a;

    if (a == 0.0) {
      std::cout << "Error: encountered a zero area in the file " << mesh_boundary_area << "\n";
      return 1;
    }

    sum_area += a;
    area_inv[i] = 1.0/a;
    area_id[i] = id;

  }

  area_inv[num_of_columns-2] = 1.0/sum_area;
  area_inv[num_of_columns-1] = 1.0/sum_area;



  std::ifstream ifs (induced_charge.c_str()); 

  std::ofstream ofs_conv ;
  std::ofstream ofs_conv_per_area;
  if (output_converted) {
    ofs_conv.open(mean_converted_induced_charge.c_str());
    ofs_conv_per_area.open(mean_converted_induced_charge_per_area.c_str());
  }

  {
    // Ignore the first line.
    std::string dummyLine;
    getline(ifs, dummyLine);

    if (output_converted) {
      ofs_conv << dummyLine << "\n";
      ofs_conv_per_area << dummyLine << "\n";
    }
  }


  std::vector<double> sum_values(num_of_columns - 1, 0);
  std::vector<double> sum_sqr_values(num_of_columns - 1, 0);  


  bool start_sum = false;
  if (start_at_timestep==-1 && start_at_line==-1) {
    start_sum = true;
  }

  unsigned num_of_values = 0;
  int read_counter = 0;  
  while (!ifs.eof()) {
    int tmp;

    ifs >> tmp;

    if (ifs.eof()) break; // don't repeat the last line

    if (output_converted) {
      ofs_conv << tmp * conversion_coef_time;
      ofs_conv_per_area << tmp * conversion_coef_time;
    }

    
    if (start_at_timestep!= -1)
      if (!start_sum && tmp >= start_at_timestep ) start_sum = true;

    if (start_at_line!= -1)
      if (!start_sum && read_counter >= start_at_line ) start_sum = true;

    for (unsigned int i = 0; i < sum_values.size(); ++i) {

      double tmp;
      ifs >> tmp;

      if (output_converted) {
        ofs_conv << " " << tmp * conversion_coef_charge;
        ofs_conv_per_area << " " << tmp * conversion_coef_charge* area_inv[i];
      }


      if (start_sum) {
        if (i == 0) ++num_of_values;   
        sum_values[i] += tmp;
      }
    }

    ofs_conv << "\n";
    ofs_conv_per_area << "\n";
    
    ++read_counter;
  }

  if (output_converted) {
    ofs_conv.close();
    ofs_conv_per_area.close();
  }
  
  if (!start_sum) {
    std::cout << "Error: start_sum has not meet the required condition\n";
    return 1;
  }

  // output file
  std::ofstream ofs (mean_induced_charge.c_str()); 
  std::ofstream ofs_area (mean_induced_charge_per_area.c_str());

  // Differential capacitance
  for (unsigned int i = 0; i < sum_values.size(); ++i) {
    double q_m = (sum_values[i]/num_of_values);
    ofs << conversion_coef_charge * q_m << " ";
    ofs_area << conversion_coef_charge*q_m*area_inv[i]<< " ";
  } 
  ofs << "\n";



}
