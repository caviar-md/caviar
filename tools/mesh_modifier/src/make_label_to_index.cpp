
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

//#include <iostream>
//#include <algorithm>    // std::sort, std::count

namespace mesh_modifier {

  void Mesh_modifier::make_label_to_index (const std::vector<Universal_dataset_number_2411> & u
                           ,std::vector<unsigned> & v) {

    unsigned max = 0;
    for (unsigned i = 0; i < u.size(); ++i) 
      if (max < u[i].record1[0]) 
        max = u[i].record1[0];
    
    v.resize(max + 1, 0); // XXX
    for (unsigned i = 0; i < u.size(); ++i) 
      v[u[i].record1[0]] = i;      
   
  }
  
  void Mesh_modifier::make_label_to_index (const std::vector<Universal_dataset_number_2412> & u
                           ,std::vector<unsigned> & v) {
    unsigned max = 0;
    for (unsigned i = 0; i < u.size(); ++i) 
      if (max < u[i].record1[0]) 
        max = u[i].record1[0];
    
    v.resize(max + 1, 0); // XXX
    for (unsigned i = 0; i < u.size(); ++i) 
      v[u[i].record1[0]] = i;            
      
  }

  
}  
