
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
// 
// this program calculates differential capacitance
// read the comments and you will know what to do.
int main () {

  // ----------------------------------------------------------
  // User variables
  // ----------------------------------------------------------

  // electrode index list. (optional).
  // since the electrode data may be divided into many parts, the user can set
  // and define the index of electrodes that belongs to one electrode.
  // for example
  // t 0 1 2 3 4 5 6 s a
  // can be divided into two electrodes
  // x 0 1 1 1 1 1 0 0 0  // electrode 1:
  // x 0 0 0 0 0 0 2 0 0  // electrode 2:
  const std::vector<std::vector<int>> electrodes_index {{0, 1, 1, 1, 1, 1, 0, 0, 0}, 
                                                        {0, 0, 0, 0, 0, 0, 1, 0, 0}};


  // If true, it will print the input file containing induced charges.
  // This is a good practice to check and compare it with the inputted file
  // to see if it is read correctly.
  const bool print_input_file = true;

  // this is a conversion coefficient optional for user. Since the simulation
  // units is usually different from report or paper units, one can convert
  // it using this variable.
  const double conversion_coef = 1.0; 


  // BOLTZMAN constant times Temperature
  const double KbT = 1.0;  


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
  const std::string diff_capacitance = "o_dc"; 

  // output file names for differential capacitance. 
  const std::string diff_capacitance_per_area = "o_dc_per_area"; 


  // output file names for differential capacitance for electrodes. 
  const std::string diff_capacitance_elec = "o_dc_elec"; 

  // output file names for differential capacitance  for electrodes. 
  const std::string diff_capacitance_per_area_elec = "o_dc_per_area_elec"; 

  // ----------------------------------------------------------
  // Calculations starts from here. 
  // ----------------------------------------------------------






  std::ifstream ifs_area (mesh_boundary_area.c_str()); 

  // importing area of the boundaries by id
  std::vector<int> area_id (num_of_columns+2,0);
  std::vector<double> area (num_of_columns+2,0);


  // electrodes part (optional)
  auto electrodes_index_size = electrodes_index.size();
  std::vector<double> area_electrodes (electrodes_index_size, 0.0);

  double sum_area = 0.0;
  for (unsigned int i =0 ; i < num_of_columns; ++i) {
    int id;
    double a;
    ifs_area >> id;    
    ifs_area >> a;
    sum_area += a;
    area[i] = a;
    area_id[i] = id;
    if (a == 0.0) {
      std::cout << "Error: encountered a zero area in the file " << mesh_boundary_area << "\n";
      return 1;
    }

    // electrodes part (optional)
    for (unsigned int j =0 ; j < electrodes_index_size; ++j) {
      area_electrodes[j] += a * electrodes_index[j][i];
    }
    
  }

  area[num_of_columns-2] = sum_area;
  area[num_of_columns-1] = sum_area;



  std::ifstream ifs (induced_charge.c_str()); 

  {
    // Ignore the first line.
    std::string dummyLine;
    getline(ifs, dummyLine);
  }

  std::vector<double> sum_values(num_of_columns - 1, 0);
  std::vector<double> sum_sqr_values(num_of_columns - 1, 0);  

  // electrodes part (optional)
  std::vector<double> sum_values_electrodes(electrodes_index_size, 0);
  std::vector<double> sum_sqr_values_electrodes(electrodes_index_size, 0);  

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

    if (print_input_file) std::cout << tmp << " ";

    if (start_at_timestep!= -1)
      if (!start_sum && tmp >= start_at_timestep ) start_sum = true;

    if (start_at_line!= -1)
      if (!start_sum && read_counter >= start_at_line ) start_sum = true;

    // electrodes part (optional)
    std::vector<double> sum_sqr_values_electrodes_tmp(electrodes_index_size, 0.0);

    for (unsigned int i = 0; i < sum_values.size(); ++i) {
        double tmp;
        ifs >> tmp;

        if (print_input_file)  std::cout <<tmp << " ";

        if (start_sum) {

          // do it once per line
          if (i == 0) ++num_of_values;   

          sum_values[i] += tmp;
          sum_sqr_values[i] += tmp*tmp;

          // electrodes part (optional)
          for (unsigned int j =0 ; j < electrodes_index_size; ++j) {
            sum_values_electrodes[j] += tmp * electrodes_index[j][i];
            sum_sqr_values_electrodes_tmp[j] += tmp * electrodes_index[j][i];
          }

        }
    }

    // electrodes part (optional)
    for (unsigned int j =0 ; j < electrodes_index_size; ++j) {
      sum_sqr_values_electrodes[j] += sum_sqr_values_electrodes_tmp[j]*sum_sqr_values_electrodes_tmp[j];
    }

    if (print_input_file) std::cout << '\n';
    
    ++read_counter;
  }
  
  if (!start_sum) {
    std::cout << "Error: start_sum has not meet the required condition\n";
    return 1;
  }

  // output file
  std::ofstream ofs (diff_capacitance.c_str()); 
  std::ofstream ofs_area (diff_capacitance_per_area.c_str());

  // Differential capacitance
  for (unsigned int i = 0; i < sum_values.size(); ++i) {
    //std::cout << sum_values[i] << " ";
    double q_m = (sum_values[i]/num_of_values);
    double q_m_sq = q_m*q_m;
    double q_sq_m =  (sum_sqr_values[i]/num_of_values);
    ofs << conversion_coef*(q_sq_m - q_m_sq)/(KbT) << " ";
    ofs_area << conversion_coef*(q_sq_m - q_m_sq)/(KbT* area[i]) << " ";
  } 
  ofs << "\n";
  ofs_area << "\n";


  // electrodes part (optional)
  std::ofstream ofs_elec (diff_capacitance_elec.c_str()); 
  std::ofstream ofs_area_elec (diff_capacitance_per_area_elec.c_str());
  for (unsigned int i =0 ; i < electrodes_index_size; ++i) {
    double q_m = (sum_values_electrodes[i]/num_of_values);
    double q_m_sq = q_m*q_m;
    double q_sq_m =  (sum_sqr_values_electrodes[i]/num_of_values);
    ofs_elec << conversion_coef*(q_sq_m - q_m_sq)/(KbT) << " ";
    ofs_area_elec << conversion_coef*(q_sq_m - q_m_sq)/(KbT* area_electrodes[i]) << " ";
  }

  ofs_elec << "\n";
  ofs_area_elec << "\n";

}
