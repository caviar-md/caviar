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


#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <map>

#include "utility.h"

//=========================================================================
//
// 1- area_counter: count number of atoms in areas at the initial and final snapshot of simulation.
// the result for atom type_x is : timestep, time, area_1_counter, area_2_counter
// each atom type results in a new file.
// there's also a file containing all of atoms.
//
// 2- transition event: record atom moving from one area to another. It may have transition to a new area and going back.
// the result for atom type_x is:
//
// 3- transition event counter: count all transition events in a matrix. It may have transition and going back.
// the result for atom type_x is the area matrix.
//
// transition counter: count complete transition from initial to final snapshot.
// the result for atom type_x is the area matrix.
//
//=========================================================================
//
// WARNING: This code may have undefined behavior for overlapping areas
//
//=========================================================================
bool debug = true;
bool debug_read = true;


// AreaMatrix example
//    | 1 | 2 | 3 |
// ----------------------
// 1  |   |   |   |
// 2  |   |   |   |
// 3  |   |   |   |
// 
// AreaMatrix is a square matrix describing transition parameters between areas
// Area 0 is defined for undefined area
class AreaMatrix
{
  public:
  AreaMatrix(){};
  
  std::vector<std::vector<int>> a;// first index: initial area. Second index: final area;
  void init(int num_areas)
  {
    a.resize(num_areas);
    for (auto && i : a) {i.resize(num_areas);}
  };
  void add_transition(int area1, int area2)
  {
    if (area1 != area2)
      a[area1][area2]++;
  }
  void clear()
  {
    a.clear();
  }
  void export_to_file(std::string filename)
  {
    int a_size = a.size();
    std::ofstream ofs;
    ofs.open(filename.c_str());
    for (int i = 0; i < a_size;++i)
    {
      for (int j = 0; j < a_size;++j)
      {
        if (j == a_size-1)
          ofs << a[i][j] << "\n";
        else
          ofs << a[i][j] << " ";
      }
    }
    ofs.close();
  }
  
  void export_to_file_csv(std::string filename)
  {
    int a_size = a.size();
    std::ofstream ofs;
    ofs.open(filename.c_str());
    for (int i = 0; i < a_size;++i)
    {
      for (int j = 0; j < a_size;++j)
      {
        if (j == a_size-1)
          ofs << a[i][j] << "\n";
        else
          ofs << a[i][j] << ",";
      }
    }
    ofs.close();
  }  
};


int get_atom_id (std::vector<int> &atom_type_to_id, int type)
{
  if (type + 1> atom_type_to_id.size()) return -1;
  return atom_type_to_id[type];
}
//=======================================================================
//=======================================================================

int main (int argc, char **argv) 
{

  //=======================================================================
  if (debug) std::cout << "---------------- Definitions -----------------\n";
  //=======================================================================
  
  std::vector < std::array <double,2> > areas1D;// intervals in 'positions_axis' direction. 1 Dimensional areas.
  areas1D.push_back({0,0}); //WARNING: Don't delete this. Index 0 is reserved for undefined area.
  
  //=======================================================================
  if (debug) std::cout << "---------------- Set Parameters -----------------\n";
  //=======================================================================
    
  std::string input_xyz_file_name = "o_xyz_to_resume.xyz";
  //std::string input_xyz_file_name = "o_xyz--.xyz";
  int xyz_steps = 10000; // each snapshot describes this number of time steps
  double time_converter = 0.000006; // time step to time unit coefficient, eg., (ps, ns, fs)
  char areas1D_axis = 'x'; // it can be 'x','y' or 'z'  
  std::vector<int> atom_types {0,1,2,3};  // the atom_types to be monitored

  //=======================================================================  
  if (debug) std::cout << "---------------- Set area intervals to be monitored -----------------\n";
  //=======================================================================

  // intervals in 'positions_axis' direction.
  // WARNING: This code may have undefined behavior for overlapping areas.
  
  areas1D.push_back({-100,-10}); 
  areas1D.push_back({-10,10}); 
  areas1D.push_back({10,100});

  //=======================================================================
  if (debug) std::cout << "---------------- Print Parameters -----------------\n";
  //=======================================================================
    
  std::cout << "input_xyz_file_name : " << input_xyz_file_name << "\n";
  std::cout << "xyz_steps : " << xyz_steps << "\n";
  std::cout << "time_converter : " << time_converter << "\n";
  std::cout << "areas1D_axis : " << areas1D_axis << "\n";
  std::cout << "the atom_types to be monitored: ";
  
  for (auto i : atom_types) std::cout << i << " , ";
  std::cout << "\n";
  
  std::cout << "Defined intervals in '" << areas1D_axis << "' axis:";
  for (int i = 0; i < areas1D.size(); ++i) {
    std::cout << " [ " << areas1D[i][0] <<" : " <<areas1D[i][1] << " ] " << ","; 
  }
  std::cout << std::endl;
  
  //=======================================================================  
  if (debug) std::cout << "---------------- Open Files-----------------\n";
  //=======================================================================

  if (!file_exists_0(input_xyz_file_name))
  {
    std::cout << "Error : '" <<input_xyz_file_name << "' file does not exist." << std::endl;
    return 1;
  }
  std::ifstream ifs (input_xyz_file_name.c_str());


  std::vector<std::shared_ptr<std::ofstream>> ofs_counter;
  for (auto &&t : atom_types)
  {
    std::shared_ptr<std::ofstream> ofs = std::make_shared<std::ofstream>();
    std::string ofs_name = "o_c_t" + std::to_string(t);
    ofs->open(ofs_name.c_str());
    ofs_counter.push_back(std::move(ofs));    
  }
  
  //=======================================================================  
  std::cout << "---------------- Initialize Variables-----------------\n";    
  //=======================================================================

  int num_areas = areas1D.size();
  int num_atom_types = atom_types.size();    
  int frame_counter = 0;    
  
  std::vector<std::vector <int> > num_atom_in_area; // first index: 'area index'. Second index : 'Atom Type id' 
  num_atom_in_area.resize(num_areas);   
  for (auto && i : num_atom_in_area)
    i.resize(num_atom_types, 0);
 
  
  /* 
  // Memory efficient 'atom_type_to_id'
  std::map<int,int> atom_type_to_id;
  //std::map<int,int> atom_id_to_type; ~ is the same as atom_types[i]
  for (int i = 0; i < num_atom_types; ++i)
  {
    type_to_id[atom_types[i]] = i;
  }
  */ 
  
  // Process Efficient 'atom_type_to_id'
  std::vector<int> atom_type_to_id(100, -1);  
  for (int i = 0; i < num_atom_types; ++i)
  {
    if (atom_type_to_id.size() < atom_types[i]+1) 
      atom_type_to_id.resize(atom_types[i]+1,-1);
    atom_type_to_id[atom_types[i]] = i;
  }

  std::vector<AreaMatrix> transition_counter      (num_atom_types);
  for (auto && c : transition_counter) c.init(num_areas + 1);
  
  std::vector<AreaMatrix> transition_event_counter(num_atom_types);
  for (auto && c : transition_event_counter) c.init(num_areas + 1);

  //std::vector<AreaMatrix> transition_event_counter(num_atom_types);

  std::vector<int> atom_type_list_in_xyz_file;      // record atoms types in xyz order
  std::vector<int> atom_area_id;      // record each atom is in which area in the current timestep; index: atom order index in xyz file.
  std::vector<int> atom_area_id_init; // record each atom is in which area at the initial timestep: index: atom order index in xyz file.
  bool initial_area = true; // used for finding first timestep
    
  //=======================================================================  
  std::cout << "---------------- File reading loop-----------------\n";
  //=======================================================================    
  
  bool first_loop = true;
  int num_of_atoms = -1;
  while (!ifs.eof()) 
  {

    //clearing container
    for (auto && x : num_atom_in_area) 
      for (auto && x2 : x)
        x2 = 0;

      
    int num_of_atoms_tmp= 0;
    
    ifs >> num_of_atoms_tmp;
    if (!first_loop)
    {
      if(num_of_atoms_tmp > 0 && num_of_atoms_tmp!=num_of_atoms)
      {
        std::cout << "Error: number of atoms changed in xyz file. Not supported. (" << num_of_atoms_tmp << " != " << num_of_atoms << ")"<< std::endl;
        return 0;
      }
    }
    num_of_atoms = num_of_atoms_tmp;
    if (atom_area_id.size() < num_of_atoms) 
    {
      atom_area_id.resize(num_of_atoms,0);
      atom_area_id_init.resize(num_of_atoms,0);
      atom_type_list_in_xyz_file.resize(num_of_atoms,-1);
    }

    if (debug_read) std::cout << "Atom " <<  num_of_atoms << "\n";
    

    if (ifs.eof()) break; // don't repeat the last line

    // Ignore the second line, 'Atom'
    std::string dummyLine;
    std::getline(ifs, dummyLine);
    if (debug_read) std::cout << "dummyLine 1: " << dummyLine << "\n";
    std::getline(ifs, dummyLine);
    if (debug_read) std::cout << "dummyLine 2: " << dummyLine << "\n";

    for (int j = 0; j < num_of_atoms; ++j) 
    {
      //===================
      // reading atom data
      //===================

      std::string line;
      std::getline(ifs, line);

      if (debug_read) std::cout << j << ") " << line << "\n";
            
      std::stringstream ss(line);  
      int type;
      ss >> type;
      atom_type_list_in_xyz_file[j]= type;      
      if (!std::count(atom_types.begin(),atom_types.end(),type)) continue;
      double x,y,z;
      ss >> x >> y >> z;


      if (debug_read) std::cout << frame_counter << " : (" << j <<" of " << num_of_atoms << ") ["<< x << "," << y << "," << z << "]\n";
      
      auto id = get_atom_id( atom_type_to_id, type);
      if (id == -1) continue; // ignore atom

      //==========================
      // processing atom location
      //==========================
      
      double tmp_coordinate = 0;
      switch(areas1D_axis)
      {
        default  :
        case 'x' : tmp_coordinate = x; break;
        case 'y' : tmp_coordinate = y; break;
        case 'z' : tmp_coordinate = z; break;
      }

      bool area_found = false;      
      for (int area_index = 1; area_index < areas1D.size(); ++area_index) 
      {

        if (areas1D[area_index][0] < tmp_coordinate && tmp_coordinate < areas1D[area_index][1]) 
        {

          //if (debug_read) std::cout << "k: " << k << " , type : " << type << " type_to_id : " 
          //                          << type_to_id[type] <<  " x_number.size " << x_number.size()
          //                          << " x_number[k].size " << x_number[k].size() << "\n"; 
          if (debug_read) std::cout << "type " << type  << " in area_1D " << area_index << "\n";

          
          num_atom_in_area[area_index][id] ++;

          if (initial_area)
          {
            atom_area_id_init[j] = area_index;
          }          
          
          //if (atom_area_id[j] != area_index)  // this will be checke in 'addd_transition) fyunction
          {
            transition_event_counter[id].add_transition (atom_area_id[j], area_index);
            atom_area_id[j] = area_index;
          }

          area_found = true;
          break;
        }                

      }

      if (!area_found) 
      {

        atom_area_id[j] = 0;
        transition_event_counter[id].add_transition(atom_area_id[j], 0);// transition to undefined area
      }

    }  
    initial_area = false;

    //-------------------------
    //writing number of atom types in areas into file
    //-------------------------
    double time = time_converter * frame_counter * xyz_steps;
    
    for (int j = 0; j < num_atom_types; ++j) 
    {

      *ofs_counter[j] << frame_counter << " " << time;

      for (int k = 0; k < num_areas; ++k)
      {
      
        *ofs_counter[j] << " " << num_atom_in_area[k][j];
      }

      *ofs_counter[j] << "\n";

    }    

    ++frame_counter;
    first_loop = false;
  }

  for (auto && ofs: ofs_counter)
  {

    ofs->close();

  }

  //=======================================================================  
  std::cout << "---------------- Writing Area Matrix into files-----------------\n";
  //=======================================================================    
  
  for (int i = 0; i<atom_area_id.size(); ++i)
  {
    auto id = get_atom_id( atom_type_to_id, atom_type_list_in_xyz_file[i]);
    if (id == -1) continue; // ignore atom
    //std::cout << i << " , id: " << id << " init_area: " << atom_area_id_init[i] << " , final_area: " << atom_area_id[i] << std::endl;
    transition_counter[id].add_transition (atom_area_id_init[i], atom_area_id[i]);    
  }
  
  for (int i = 0 ; i < transition_counter.size(); ++i)
  {
    std::string filename = "o_tc_"+std::to_string(atom_types[i]);
    transition_counter[i].export_to_file(filename);
  }
  
  
  for (int i = 0 ; i < transition_event_counter.size(); ++i)
  {
    std::string filename = "o_tec_"+std::to_string(atom_types[i]);
    transition_event_counter[i].export_to_file(filename);
  }  
  
  std::cout << "number of xyz frames: "<< frame_counter << std::endl;

}

//std::cout << "argc: " << argc << "\nargv:";
/*
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
*/
