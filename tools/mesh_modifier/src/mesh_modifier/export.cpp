
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

namespace mesh_modifier {

  void Mesh_modifier::export_ (const std::string & filename, bool unsupported) {

    unsigned uci = unv_container.size() - 1;
    if (uci<0) {
      std::cout << "error: there's no unv_container\n";
      return;
    }
    
    std::ofstream ofs;
    ofs.open (filename.c_str());
        
    if (unsupported) {
      if (unv_container[uci].udn_unsupported.size()>0)
        export_udn_unsupported (ofs);
    }
    
    export_udn_2411 (ofs);
    export_udn_2412 (ofs);
    export_udn_2467 (ofs);
    
    ofs.close ();    
  }
  
  void Mesh_modifier::export_udn_unsupported (std::ofstream & ofs) {
    unsigned uci = unv_container.size() - 1;  
//    ofs << std::right << std::setw(6) << "-1" << "\n";
    for (unsigned i = 0; i < unv_container[uci].udn_unsupported.size(); ++i) {
      ofs << std::right << std::setw(6) << "-1" << "\n";          
      for (unsigned j = 0; j < unv_container[uci].udn_unsupported[i].records.size(); ++j) {
        ofs << unv_container[uci].udn_unsupported[i].records[j];
      }
      ofs << std::right << std::setw(6) << "-1" << "\n";      
    }
    
  }
  
    
  void Mesh_modifier::export_udn_2411 (std::ofstream & ofs) {
    unsigned uci = unv_container.size() - 1;  
    ofs << std::right << std::setw(6) << "-1" << "\n";               
    ofs << std::right << std::setw(6) << "2411" << "\n";                   
    for (unsigned i = 0; i < unv_container[uci].udn_2411.size(); ++i) {
      ofs << std::right 
          << std::setw(10) << unv_container[uci].udn_2411[i].record1[0]
          << std::setw(10) << unv_container[uci].udn_2411[i].record1[1]
          << std::setw(10) << unv_container[uci].udn_2411[i].record1[2]
          << std::setw(10) << unv_container[uci].udn_2411[i].record1[3] << "\n";      
      

      ofs << std::scientific << std::setprecision(16) << std::uppercase
          << std::setw(25) << unv_container[uci].udn_2411[i].record2[0] 
          << std::setw(25) << unv_container[uci].udn_2411[i].record2[1]
          << std::setw(25) << unv_container[uci].udn_2411[i].record2[2] << "\n"; 
      ofs << std::fixed;
           
    }
    ofs << std::right << std::setw(6) << "-1" << "\n";          
    
  }
  

  void Mesh_modifier::export_udn_2412 (std::ofstream & ofs) {
    unsigned uci = unv_container.size() - 1;  
    ofs << std::right << std::setw(6) << "-1" << "\n";               
    ofs << std::right << std::setw(6) << "2412" << "\n";
                       
    for (unsigned i = 0; i < unv_container[uci].udn_2412.size(); ++i) {
    
      ofs << std::right 
          << std::setw(10) << unv_container[uci].udn_2412[i].record1[0]
          << std::setw(10) << unv_container[uci].udn_2412[i].record1[1]
          << std::setw(10) << unv_container[uci].udn_2412[i].record1[2]
          << std::setw(10) << unv_container[uci].udn_2412[i].record1[3]
          << std::setw(10) << unv_container[uci].udn_2412[i].record1[4]
          << std::setw(10) << unv_container[uci].udn_2412[i].record1[5] << "\n";  
              
      int FE_Id = unv_container[uci].udn_2412[i].record1[1];
      bool beam_type = (FE_Id==11 || FE_Id==21 || FE_Id==22 || FE_Id==23 || FE_Id==24);
      unsigned num_of_elements = unv_container[uci].udn_2412[i].record1[5];
      
      if (beam_type) {
      
        ofs << std::setw(10) << unv_container[uci].udn_2412[i].record2[0]
            << std::setw(10) << unv_container[uci].udn_2412[i].record2[1]
            << std::setw(10) << unv_container[uci].udn_2412[i].record2[2] << "\n";
            
        for (unsigned int j = 0; j < num_of_elements; ++j)     
          ofs << std::setw(10) << unv_container[uci].udn_2412[i].record3[j];
        ofs << "\n";              
        
      } else {
      
        for (unsigned int j = 0; j < num_of_elements; ++j)     
          ofs << std::setw(10) << unv_container[uci].udn_2412[i].record2[j];      
        ofs << "\n";                    
        
      }
    }
    ofs << std::right << std::setw(6) << "-1" << "\n";  
  }                  
   
                

  void Mesh_modifier::export_udn_2467 (std::ofstream & ofs) {
    unsigned uci = unv_container.size() - 1;
    if (unv_container[uci].udn_2467.size() == 0)
      return;
    ofs << std::right << std::setw(6) << "-1" << "\n";               
    ofs << std::right << std::setw(6) << "2467" << "\n";                   
    for (unsigned i = 0; i < unv_container[uci].udn_2467.size(); ++i) {  
    

      ofs << std::right 
          << std::setw(10) << unv_container[uci].udn_2467[i].record1[0]
          << std::setw(10) << unv_container[uci].udn_2467[i].record1[1]
          << std::setw(10) << unv_container[uci].udn_2467[i].record1[2]
          << std::setw(10) << unv_container[uci].udn_2467[i].record1[3]
          << std::setw(10) << unv_container[uci].udn_2467[i].record1[4]
          << std::setw(10) << unv_container[uci].udn_2467[i].record1[5]
          << std::setw(10) << unv_container[uci].udn_2467[i].record1[6]
          << std::setw(10) << unv_container[uci].udn_2467[i].record1[7] << "\n";
                                                                                                
     //group name
      ofs << unv_container[uci].udn_2467[i].record2 << "\n";
      
      unsigned num_of_elements = unv_container[uci].udn_2467[i].record1[7];
      
      bool print_eol = false;
      for (unsigned j = 0; j < num_of_elements*4; ++j) {
        ofs << std::setw(10) << unv_container[uci].udn_2467[i].record3[j];
        print_eol = true;              
        if (((j+1)%8)==0) {
          ofs << "\n";
          print_eol = false;                  
        }
      }
      if (print_eol) ofs << "\n";

    }   
    
    ofs << std::right << std::setw(6) << "-1" << "\n";  
  } 
}


