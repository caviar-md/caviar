
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

#include <iostream>

namespace mesh_modifier {

  void Mesh_modifier::merge_all (bool unsupported) {
  
    if (unv_container.size()<2) {
      std::cout << "error: there has to be at least two unv_container to merge.\n";
      return;
    }
    
    unv_container.push_back(Unv_container());
    
    unsigned current_container = unv_container.size() - 1;
    auto &c1  = unv_container[current_container];
    
    
    for (unsigned i = 0; i < current_container; ++i) {
      add_to_unv_container (c1, unv_container[i]);
    }
    
    if (unsupported) {
      std::cout << "warning: only the unsupported udn sections" 
                << "of the first imported file will be added.\n";
                
      c1.udn_unsupported.insert(std::end(c1.udn_unsupported), 
                         std::begin(unv_container[0].udn_unsupported),
                         std::end(unv_container[0].udn_unsupported));
    }
    
    fix_udn_2411_element_labels ();
    fix_udn_2412_element_labels ();    
    fix_udn_2412_node_labels ();
    fix_udn_2467_element_labels ();        
    fix_udn_2467_entity_tags ();        
  }
  
  void Mesh_modifier::add_to_unv_container (Unv_container &c1, Unv_container &c2) {
    c1.udn_2411.insert(std::end(c1.udn_2411), std::begin(c2.udn_2411), std::end(c2.udn_2411));
    c1.udn_2412.insert(std::end(c1.udn_2412), std::begin(c2.udn_2412), std::end(c2.udn_2412)); 
    c1.udn_2467.insert(std::end(c1.udn_2467), std::begin(c2.udn_2467), std::end(c2.udn_2467));    
  }

  void Mesh_modifier::fix_udn_2411_element_labels () {  
    unsigned uci = unv_container.size() - 1;    
    
    unsigned k = 0;
    unsigned previous_size = 0;
    
    for (unsigned i = 0; i < uci; ++i) {
      if (i != 0 ) previous_size += unv_container[i-1].udn_2411.size();

      for (unsigned j = 0; j < unv_container[i].udn_2411.size(); ++j) {
        unv_container[uci].udn_2411[k].record1[0] += previous_size;
        ++k;

      }
    }    
  }
  
   void Mesh_modifier::fix_udn_2412_element_labels () {  
    unsigned uci = unv_container.size() - 1;    
    
    unsigned k = 0;
    unsigned previous_size = 0;
    
    for (unsigned i = 0; i < uci; ++i) {
      if (i != 0 ) previous_size += unv_container[i-1].udn_2412.size();

      for (unsigned j = 0; j < unv_container[i].udn_2412.size(); ++j) {
        unv_container[uci].udn_2412[k].record1[0] += previous_size;
        ++k;
      }
    }    
  } 

  void Mesh_modifier::fix_udn_2412_node_labels () {  
    unsigned uci = unv_container.size() - 1;    
    
    unsigned k = 0;
    unsigned previous_size = 0;
    
    for (unsigned i = 0; i < uci; ++i) {
      if (i != 0 ) previous_size += unv_container[i-1].udn_2411.size(); //size of point labels//

      for (unsigned j = 0; j < unv_container[i].udn_2412.size(); ++j) {
        //unv_container[uci].udn_2412[k].record1[0] += previous_size;
        int FE_Id = unv_container[uci].udn_2412[k].record1[1];
        bool beam_type = (FE_Id==11 || FE_Id==21 || FE_Id==22 || FE_Id==23 || FE_Id==24);
        unsigned num_of_elements = unv_container[uci].udn_2412[k].record1[5];
                            
        if (beam_type) { //beam elements
          for (unsigned int m = 0; m < num_of_elements; ++m)
            unv_container[uci].udn_2412[k].record3[m] += previous_size;
        } else {
          for (unsigned int m = 0; m < num_of_elements; ++m)
            unv_container[uci].udn_2412[k].record2[m] += previous_size;        
        }
        
        ++k;
      }
    }    
  } 

  void Mesh_modifier::fix_udn_2467_element_labels () {  
    unsigned uci = unv_container.size() - 1;    
    
    unsigned k = 0;
    unsigned previous_size = 0;
    
    for (unsigned i = 0; i < uci; ++i) {
      if (i != 0 ) previous_size += unv_container[i-1].udn_2467.size();

      for (unsigned j = 0; j < unv_container[i].udn_2467.size(); ++j) {
        unv_container[uci].udn_2467[k].record1[0] += previous_size;
        ++k;

      }
    }    
  }
    
  void Mesh_modifier::fix_udn_2467_entity_tags () {  
    unsigned uci = unv_container.size() - 1;    
    

    unsigned s = 0;
    unsigned previous_size = 0;    
    for (unsigned i = 0; i < uci; ++i) {
      if (i != 0 ) previous_size += unv_container[i-1].udn_2412.size();

      for (unsigned j = 0; j < unv_container[i].udn_2467.size(); ++j) {

        unsigned num_of_entities = unv_container[i].udn_2467[j].record1[7];

        for (unsigned k = 0; k < num_of_entities; ++k) {
          unsigned m = (4*k) + 1;
          unv_container[uci].udn_2467[s].record3[m] += previous_size;

        }
        ++s;
      }
    }    
  }
}







