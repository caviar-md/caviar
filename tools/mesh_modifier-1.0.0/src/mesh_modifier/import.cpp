
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

namespace mesh_modifier {


   
  
  void Mesh_modifier::import (const std::string & filename, bool unsupported) {
    unv_container.push_back(Unv_container());
    
    std::ifstream ifs;
    ifs.open (filename.c_str());
    bool ignore_section = false;
    bool in_section = false;
    int udn_code = 0;
    while (!ifs.eof()) {
      if (!in_section) {
        std::string c;
        ifs >> c;
        if (c == "#") 
          ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if (c == "-1") {
//        std::cout << "\n";
          in_section = true;
          ifs >> udn_code;
          std::cout << "udn code : " << udn_code << "\n";
        }
      } else {
        if (udn_code == 2411) import_udn_2411 (ifs);
        else if (udn_code == 2412) import_udn_2412 (ifs);
        else if (udn_code == 2467) import_udn_2467 (ifs);
        else {

          if (unsupported)
            import_udn_unsupported (ifs, udn_code);
          else
            import_udn_ignore (ifs, udn_code);
        }    
        in_section = false;
      }
    }
    ifs.close ();
  }

  void Mesh_modifier::import_udn_ignore (std::ifstream & ifs, int udn_code) {
    std::cout << "Warning: Unsupported udn " << udn_code 
              << ". Ignoring the section.\n";  
    while (true) {
     std::string c;
     ifs >> c;
     if (c == "-1") return;
    }  
  }
  
  void Mesh_modifier::import_udn_unsupported (std::ifstream & ifs, int udn_code) {
    std::cout << "Warning: Unsupported udn " << udn_code 
              << ". Importing it as a string.\n";  
    unsigned uci = unv_container.size() - 1;
    unv_container[uci].udn_unsupported.push_back (Universal_dataset_number_unsupported());
    unsigned udi = unv_container[uci].udn_unsupported.size() - 1;
    std::stringstream a_stream;
    a_stream << std::right << std::setw(6) << udn_code << "\n";
    unv_container[uci].udn_unsupported[udi].records.push_back (a_stream.str());    
    while (true)
    {
      std::string line;
      std::getline(ifs, line);
      std::vector<std::string> tokens;
      std::istringstream iss (line);
      std::copy (std::istream_iterator<std::string>(iss)
                ,std::istream_iterator<std::string>()
                ,std::back_inserter(tokens));
      if (tokens.size()>0) {
        if (tokens[0]=="-1")
          break;       
      line.append("\n");  
      unv_container[uci].udn_unsupported[udi].records.push_back (line);
      }    

    } 

  }
  
  void Mesh_modifier::import_udn_2411 (std::ifstream & ifs) {
    unsigned uci = unv_container.size() - 1;   
    while (true) {
      int tmp;
      ifs >> tmp;
      if (tmp == -1) return;

      unv_container[uci].udn_2411.push_back (Universal_dataset_number_2411());
      unsigned udi = unv_container[uci].udn_2411.size() - 1;
      
      
      int dummy[3];
      ifs >> dummy[0] >> dummy[1] >> dummy[2];
      //std::cout<< tmp << " " <<  dummy[0] << " " << dummy[1] << " " << dummy[2] << "\n";
      
      unv_container[uci].udn_2411[udi].record1[0] = tmp;
      unv_container[uci].udn_2411[udi].record1[1] = dummy[0];
      unv_container[uci].udn_2411[udi].record1[2] = dummy[1];
      unv_container[uci].udn_2411[udi].record1[3] = dummy[2];            
      
      double pos[3];      
      ifs >> pos[0] >> pos[1] >> pos[2];
      //std::cout << pos[0] << " " << pos[1] << " " << pos[2] << "\n"; 
      
      unv_container[uci].udn_2411[udi].record2[0] = pos[0];
      unv_container[uci].udn_2411[udi].record2[1] = pos[1];
      unv_container[uci].udn_2411[udi].record2[2] = pos[2];

    }
     
  }

  void Mesh_modifier::import_udn_2412 (std::ifstream & ifs) {
    unsigned uci = unv_container.size() - 1;
    while (true) {
      int tmp;
      ifs >> tmp;
      if (tmp == -1) return;

      unv_container[uci].udn_2412.push_back (Universal_dataset_number_2412());
      unsigned udi = unv_container[uci].udn_2412.size() - 1;
      
      int dummy[5];
      ifs >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> dummy[4];
      
      //std::cout << tmp << " " <<  dummy[0] << " " << dummy[1] << " " << dummy[2] 
      //          << " " <<  dummy[3] << " " << dummy[4] << "\n";
      unv_container[uci].udn_2412[udi].record1.push_back( tmp);
      unv_container[uci].udn_2412[udi].record1.push_back( dummy[0]);
      unv_container[uci].udn_2412[udi].record1.push_back( dummy[1]);
      unv_container[uci].udn_2412[udi].record1.push_back( dummy[2]);       
      unv_container[uci].udn_2412[udi].record1.push_back( dummy[3]);
      unv_container[uci].udn_2412[udi].record1.push_back( dummy[4]);
      
      int FE_Id = dummy[0];
      bool beam_type = (FE_Id==11 || FE_Id==21 || FE_Id==22 || FE_Id==23 || FE_Id==24);
      unsigned num_of_elements = dummy[4];
                            
      if (beam_type) { //beam elements
        int dummy_[3];
        ifs >> dummy_[0] >> dummy_[1] >> dummy_[2];
        //std::cout <<  dummy_[0] << " " << dummy_[1] << " " << dummy_[2]  << "\n";
        unv_container[uci].udn_2412[udi].record2.push_back(dummy_[0]);
        unv_container[uci].udn_2412[udi].record2.push_back(dummy_[1]);
        unv_container[uci].udn_2412[udi].record2.push_back(dummy_[2]);
        
        int field[num_of_elements];
        for (unsigned int i = 0; i < num_of_elements; ++i)
          ifs >> field[i];
        for (unsigned int i = 0; i < num_of_elements; ++i)
          //std::cout << field[i] << " ";
          unv_container[uci].udn_2412[udi].record3.push_back(field[i]);
        //std::cout << "\n";
         
      } else {
        int field[num_of_elements];
        for (unsigned int i = 0; i < num_of_elements; ++i)
          ifs >> field[i];
        for (unsigned int i = 0; i < num_of_elements; ++i)
          //std::cout << field[i] << " ";
          unv_container[uci].udn_2412[udi].record2.push_back(field[i]);
        //std::cout << "\n";          
      }

    }
     
  }
  
  
  void Mesh_modifier::import_udn_2467 (std::ifstream & ifs) {
    unsigned uci = unv_container.size() - 1;
    while (true) {
 
      int tmp;
      ifs >> tmp;
      if (tmp == -1) return;

      unv_container[uci].udn_2467.push_back (Universal_dataset_number_2467());
      unsigned udi = unv_container[uci].udn_2467.size() - 1;
      
      int dummy[7];
      ifs >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> dummy[4]
          >> dummy[5] >> dummy[6] ;
      //std::cout << tmp << " " <<  dummy[0] << " " << dummy[1] << " " 
      //          << dummy[2]   << " " << dummy[3] << " " << dummy[4] << " " 
      //          << dummy[5]   << " " << dummy[6] << "\n";
      
      unv_container[uci].udn_2467[udi].record1.push_back(tmp);     
      for (unsigned int i = 0; i < 7; ++i)
        unv_container[uci].udn_2467[udi].record1.push_back (dummy[i]);
        
      int group_name;
      ifs >> group_name;
//      std::cout << group_name << "\n";
      unv_container[uci].udn_2467[udi].record2 = group_name;     
      
      bool print_eol = false;
      for (int i = 0; i < dummy[6]; ++i) {
        int dummy[4];
        ifs >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] ;        
//        std::cout << dummy[0] << " " << dummy[1] << " " << dummy[2] << " " 
//                  << dummy[3] ;
        for (unsigned int i = 0; i < 4; ++i)
          unv_container[uci].udn_2467[udi].record3.push_back (dummy[i]);
          
        if (print_eol) {
          //std::cout << "\n";
          print_eol = false;
        } else {
          //std::cout << " ";
          print_eol = true;      
        }
      }
      if (print_eol) std::cout << "\n";

    }
  }
}


