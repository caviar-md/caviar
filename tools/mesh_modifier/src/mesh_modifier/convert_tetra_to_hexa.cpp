
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

#include "vector.h"

#include <iostream>
#include <algorithm>


namespace mesh_modifier {

  void Mesh_modifier::convert_tetra_to_hexa () {
    unsigned uci = unv_container.size() - 1;
    if (uci<0) {
      std::cout << "error: there's no unv_container\n";
      return;
    }
    
//  sort_all_labels (); TODO later
    std::cout << "find_all_boundary_and_material_ids:" << std::endl;
    
    std::vector<unsigned> udn_2467_ids (unv_container[uci].udn_2412.size(), -1);
    extract_udn_2467_index_ids (udn_2467_ids);

    std::vector<unsigned> tetra_edges_index;
    
    
    std::cout << "make_hexa_edge_from_tetra_edge:" << std::endl;             
    make_hexa_edge_from_tetra_edge (udn_2467_ids, tetra_edges_index);
    std::cout << "make_quad_from_tiangle:" << std::endl;             
    make_quad_from_triangle (udn_2467_ids);
    std::cout << "make_hexa_from_tetra:" << std::endl;             
    make_hexa_from_tetra (udn_2467_ids);
    

    std::cout << "remove_all_tetra_edges:" << std::endl;
    remove_all_tetra_edges (tetra_edges_index);
    
    std::cout << "remove_all_triangles:" << std::endl;    
    remove_all_triangles ();
    
    std::cout << "remove_all_tetrahedrons:" << std::endl;    
    remove_all_tetrahedrons ();    
    

    
  }

  void Mesh_modifier::remove_all_tetra_edges (std::vector<unsigned> & tetra_edges) {

    unsigned uci = unv_container.size() - 1;
    auto & udn_2412 = unv_container[uci].udn_2412;/*
    for (unsigned int i = 0; i < tetra_edges.size(); ++i) {
        udn_2412.erase(udn_2412.begin() + tetra_edges[i]);    
    }*/
    for (int i = udn_2412.size() - 1 ; i > -1; --i) {
      if (udn_2412[i].record1[1] == 11)
        udn_2412.erase(udn_2412.begin() + i);
    }      
  } 
 
  
  void Mesh_modifier::remove_all_triangles () {
    unsigned uci = unv_container.size() - 1;
    auto & udn_2412 = unv_container[uci].udn_2412;
    for (int i = udn_2412.size() - 1 ; i > -1; --i) {
      if (udn_2412[i].record1[1] == 41)
        udn_2412.erase(udn_2412.begin() + i);
    }    
  }
  
  void Mesh_modifier::remove_all_tetrahedrons () {
    unsigned uci = unv_container.size() - 1;
    auto & udn_2412 = unv_container[uci].udn_2412;
    for (int i = udn_2412.size() - 1 ; i > -1; --i) {
      if (udn_2412[i].record1[1] == 111)
        udn_2412.erase(udn_2412.begin() + i);
    }
  }  

  void Mesh_modifier::extract_udn_2467_index_ids (std::vector<unsigned> & ids) {
                                   
    unsigned uci = unv_container.size() - 1;
    unsigned num_of_2467 = unv_container[uci].udn_2467.size();
    auto & udn_2467 = unv_container[uci].udn_2467;    

    std::vector<unsigned> li;
    make_label_to_index (unv_container[uci].udn_2412, li);  

    for (unsigned i = 0 ; i < num_of_2467; ++i) {
      unsigned num_of_elements = udn_2467[i].record1[7];    
      for (unsigned j = 0; j < num_of_elements; ++j) {
        unsigned m = 4*j + 1;
        unsigned label = udn_2467[i].record3[m];
        ids [li[label]] = i;
      }
    }    
    
  }



  
  void Mesh_modifier::make_hexa_edge_from_tetra_edge (std::vector<unsigned> & udn_2467_ids,
                                                    std::vector<unsigned> & tetra_edges) {
    unsigned uci = unv_container.size() - 1;
    unsigned num_of_2412 = unv_container[uci].udn_2412.size();
    auto & udn_2411 = unv_container[uci].udn_2411;
    for (unsigned i = 0 ; i < num_of_2412; ++i) {
      auto & udn_2412 = unv_container[uci].udn_2412[i];
      if (udn_2412.record1[1] == 11) {
        tetra_edges.push_back(i);
//      XXX : Spliting tetrahedrons to hexahedrons makes all the boundary lines
//            disappear. So creating them won't help. 

        /*
          tetrahedron point lables
        */        
//        auto pl0 = udn_2412.record3[0]; 
//        auto pl1 = udn_2412.record3[1];
        /*
         the index may be different from the label.
        */                
//        auto pi0 = index_of_point (pl0);
//        auto pi1 = index_of_point (pl1);        
        /*
          vertices of the tetrahedron
        */        
//        Vector<double> p0 {udn_2411[pi0].record2[0], udn_2411[pi0].record2[1], udn_2411[pi0].record2[2]};
//        Vector<double> p1 {udn_2411[pi1].record2[0], udn_2411[pi1].record2[1], udn_2411[pi1].record2[2]};        
        /*
          edges, triangles and the tetrahedron central points
        */
//        Vector<double> p01 = (p0 + p1)*0.5;        
        /*
          add the points and get its labels
        */        
        
//        auto pl01 = add_point_to_2411 (p01);     
        
        /*
          lines (beams or edges) connecting the center of triangles to vertices,
          some of them may become internal and have to be removed.
          We won't create explicit internal lines.
        */        

//        auto ll1 = add_edge_to_2412 (pl0, pl01); 
//        auto ll2 = add_edge_to_2412 (pl01, pl1);          

        /*
          add material ids to udn_2467
        */
        if (udn_2467_ids[i] != -1) {
//          add_label_to_udn_2467 (ll1, udn_2467_ids[i]); 
//          add_label_to_udn_2467 (ll2, udn_2467_ids[i]); 
          remove_label_from_udn_2467 (udn_2412.record1[0], udn_2467_ids[i]);          
        }
                                                
      }
    }
  
  }

  
  void Mesh_modifier::make_quad_from_triangle (std::vector<unsigned> & udn_2467_ids) {
    unsigned uci = unv_container.size() - 1;
    unsigned num_of_2412 = unv_container[uci].udn_2412.size();
    auto & udn_2411 = unv_container[uci].udn_2411;
    for (unsigned i = 0 ; i < num_of_2412; ++i) {
      auto & udn_2412 = unv_container[uci].udn_2412[i];
      if (udn_2412.record1[1] == 41) {


        //
        //  triangle point lables
        //        
        auto pl0 = udn_2412.record2[0]; 
        auto pl1 = udn_2412.record2[1];
        auto pl2 = udn_2412.record2[2];

        //
        // the index may be different from the label.
        //                
        auto pi0 = index_of_point (pl0);
        auto pi1 = index_of_point (pl1);
        auto pi2 = index_of_point (pl2);

        //
        // vertices of the tetrahedron
        //        
        Vector<double> p0 {udn_2411[pi0].record2[0], udn_2411[pi0].record2[1], udn_2411[pi0].record2[2]};
        Vector<double> p1 {udn_2411[pi1].record2[0], udn_2411[pi1].record2[1], udn_2411[pi1].record2[2]};
        Vector<double> p2 {udn_2411[pi2].record2[0], udn_2411[pi2].record2[1], udn_2411[pi2].record2[2]};

        //
        //  edges, triangles and the tetrahedron central points
        //
        Vector<double> p01 = (p0 + p1)*0.5;
        Vector<double> p02 = (p0 + p2)*0.5;
        Vector<double> p12 = (p1 + p2)*0.5;
                
        const double one_third = 1.0 / 3.0;
        Vector<double> p012 = (p0 + p1 + p2)*one_third;     

        std::cout << i << " " <<std::flush;
        //
        //  add the points and get its labels
        //        
        auto pl01 = add_point_to_2411 (p01);
        auto pl02 = add_point_to_2411 (p02);        
        auto pl12 = add_point_to_2411 (p12);                        
        auto pl012 = add_point_to_2411 (p012);                     

           
        //
        //  quads created from the division of the tetrahedron triangle,
        //  some of them may become internal and have to be removed.
        //  We won't create explicit internal quad.          
        //
        auto ql1 = add_quad_to_2412 (pl012, pl01, pl1, pl12); 
        auto ql2 = add_quad_to_2412 (pl012, pl02, pl0, pl01); 
        auto ql3 = add_quad_to_2412 (pl012, pl12, pl2, pl02);             

        //
        //  add material ids to udn_2467
        //
        if (udn_2467_ids[i] != -1) {
          add_label_to_udn_2467 (ql1, udn_2467_ids[i]);
          add_label_to_udn_2467 (ql2, udn_2467_ids[i]);
          add_label_to_udn_2467 (ql3, udn_2467_ids[i]);
          remove_label_from_udn_2467 (udn_2412.record1[0], udn_2467_ids[i]);          
        }                                   
            
      }
    }
  
  }

  
  void Mesh_modifier::make_hexa_from_tetra (std::vector<unsigned> & udn_2467_ids) {
    unsigned uci = unv_container.size() - 1;
    unsigned num_of_2412 = unv_container[uci].udn_2412.size();
    auto & udn_2411 = unv_container[uci].udn_2411;
    for (unsigned i = 0 ; i < num_of_2412; ++i) {
      auto & udn_2412 = unv_container[uci].udn_2412[i];
      if (udn_2412.record1[1] == 111) { 
      
        //
        //  tetrahedron point lables
        //
        
        auto pl0 = udn_2412.record2[0]; 
        auto pl1 = udn_2412.record2[1];
        auto pl2 = udn_2412.record2[2];
        auto pl3 = udn_2412.record2[3];

        //
        // the index may be different from the label.
        //                
        auto pi0 = index_of_point (pl0);
        auto pi1 = index_of_point (pl1);
        auto pi2 = index_of_point (pl2);
        auto pi3 = index_of_point (pl3);                             

        
        //
        //  vertices of the tetrahedron
        //        
        Vector<double> p0 {udn_2411[pi0].record2[0], udn_2411[pi0].record2[1], udn_2411[pi0].record2[2]};
        Vector<double> p1 {udn_2411[pi1].record2[0], udn_2411[pi1].record2[1], udn_2411[pi1].record2[2]};
        Vector<double> p2 {udn_2411[pi2].record2[0], udn_2411[pi2].record2[1], udn_2411[pi2].record2[2]};
        Vector<double> p3 {udn_2411[pi3].record2[0], udn_2411[pi3].record2[1], udn_2411[pi3].record2[2]};
        
        //
        //  edges, triangles and the tetrahedron central points
        //
        Vector<double> p01 = (p0 + p1)*0.5;
        Vector<double> p02 = (p0 + p2)*0.5;
        Vector<double> p03 = (p0 + p3)*0.5;
        Vector<double> p12 = (p1 + p2)*0.5;
        Vector<double> p13 = (p1 + p3)*0.5;
        Vector<double> p23 = (p2 + p3)*0.5;
        
        const double one_third = 1.0 / 3.0;
        Vector<double> p012 = (p0 + p1 + p2)*one_third;
        Vector<double> p013 = (p0 + p1 + p3)*one_third;
        Vector<double> p023 = (p0 + p2 + p3)*one_third;
        Vector<double> p123 = (p1 + p2 + p3)*one_third;                                               
        
        Vector<double> p0123 = (p0+ p1 + p2 + p3)*0.25;
        
        //
        //  add the points and get its labels
        //            
        auto pl01 = add_point_to_2411 (p01);
        auto pl02 = add_point_to_2411 (p02);        
        auto pl03 = add_point_to_2411 (p03);        
        auto pl12 = add_point_to_2411 (p12);
        auto pl13 = add_point_to_2411 (p13);
        auto pl23 = add_point_to_2411 (p23);   
                        
        auto pl012 = add_point_to_2411 (p012);
        auto pl013 = add_point_to_2411 (p013);        
        auto pl023 = add_point_to_2411 (p023);        
        auto pl123 = add_point_to_2411 (p123);

        auto pl0123 = add_point_to_2411 (p0123);         

        //
        //   hexahedrons created from points
        //
        
        auto hl1 = add_hexahedron_to_2412 (pl0123, pl023, pl23, pl123
                                          ,pl012 , pl02 , pl2 , pl12 ); 
        auto hl2 = add_hexahedron_to_2412 (pl0123, pl013, pl13, pl123
                                          ,pl023 , pl03 , pl3 , pl23 );
        auto hl3 = add_hexahedron_to_2412 (pl0123, pl012, pl12, pl123
                                          ,pl013 , pl01 , pl1 , pl13 );
        auto hl4 = add_hexahedron_to_2412 (pl0123, pl023, pl02, pl012
                                          ,pl013 , pl03,  pl0 , pl01); 
        //
        //  add material ids to udn_2467
        //
        if (udn_2467_ids[i] != -1) {
          add_label_to_udn_2467 (hl1, udn_2467_ids[i]);
          add_label_to_udn_2467 (hl2, udn_2467_ids[i]);
          add_label_to_udn_2467 (hl3, udn_2467_ids[i]);
          add_label_to_udn_2467 (hl4, udn_2467_ids[i]);
          remove_label_from_udn_2467 (udn_2412.record1[0], udn_2467_ids[i]);          
        }
                                               
      }
    }
  
  }
  

  void Mesh_modifier::remove_label_from_udn_2467 (unsigned element_label
                                               ,unsigned group_id) {
    unsigned uci = unv_container.size() - 1;
    auto & udn_2467 = unv_container[uci].udn_2467[group_id];
    for (unsigned int i = 0; i < udn_2467.record1[7]; ++i) {
      unsigned m = 4*i + 1;
      if (element_label == udn_2467.record3[m]) {
        udn_2467.record3.erase(udn_2467.record3.begin() + m+2);
        udn_2467.record3.erase(udn_2467.record3.begin() + m+1);
        udn_2467.record3.erase(udn_2467.record3.begin() + m);
        udn_2467.record3.erase(udn_2467.record3.begin() + m-1);        
        udn_2467.record1[7] -= 1;
        return;
      }
    }

    std::cout << "Error : element_label " << element_label 
              << " could not be found at udn_2467\n"; 
  }  

  void Mesh_modifier::add_label_to_udn_2467 (unsigned element_label
                                              ,unsigned group_id) {
    unsigned uci = unv_container.size() - 1;

    auto & udn_2467 = unv_container[uci].udn_2467[group_id];        
    udn_2467.record3.push_back(8);
    udn_2467.record3.push_back(element_label);
    udn_2467.record3.push_back(0);
    udn_2467.record3.push_back(0);
    udn_2467.record1[7] += 1;    
  }  
  
  unsigned Mesh_modifier::index_of_point (unsigned label) {
    unsigned uci = unv_container.size() - 1;
    auto & udn_2411 = unv_container[uci].udn_2411;
    for (unsigned int i = 0; i <  udn_2411.size(); ++i)
      if (udn_2411[i].record1[0] == label)
        return i; 
    std::cout << "Error: could not find a point with label." << label << "\n";
    return 0;
  }
  
  unsigned Mesh_modifier::add_point_to_2411 (Vector<double>& v) {

    unsigned uci = unv_container.size() - 1;

    auto & udn_2411 = unv_container[uci].udn_2411;
//--- new
    auto indx = udn_2411.size();
    int label = 1;
    if (indx>0) label = udn_2411[indx-1].record1[0] + 1; 

    Universal_dataset_number_2411 u;
    u.record1[0] = label;
    u.record1[1] = 1;
    u.record1[2] = 1;
    u.record1[3] = 11;
    u.record2[0] = v.x;
    u.record2[1] = v.y;
    u.record2[2] = v.z; 

    udn_2411.push_back(u);

// old
/*
    udn_2411.push_back(Universal_dataset_number_2411());

    auto indx = udn_2411.size() - 1;
    auto last_label = udn_2411[indx-1].record1[0];
    auto label = last_label + 1;

    udn_2411[indx].record1[0] = label;
    udn_2411[indx].record1[1] = 1;
    udn_2411[indx].record1[2] = 1;
    udn_2411[indx].record1[3] = 11;
    udn_2411[indx].record2[0] = v.x;
    udn_2411[indx].record2[1] = v.y;
    udn_2411[indx].record2[2] = v.z;
*/
    return label;    
  }
  
  unsigned Mesh_modifier::add_edge_to_2412 (unsigned l1, unsigned l2) {
    unsigned uci = unv_container.size() - 1;
    auto & udn_2412 = unv_container[uci].udn_2412;
    udn_2412.push_back(Universal_dataset_number_2412());
    auto indx = udn_2412.size() - 1;
    auto last_label = udn_2412[indx-1].record1[0];
    auto label = last_label + 1;    
    udn_2412[indx].record1.push_back (label);
    udn_2412[indx].record1.push_back (11);
    udn_2412[indx].record1.push_back (2);
    udn_2412[indx].record1.push_back (1);
    udn_2412[indx].record1.push_back (7);
    udn_2412[indx].record1.push_back (2);
    udn_2412[indx].record2.push_back (0);
    udn_2412[indx].record2.push_back (0);
    udn_2412[indx].record2.push_back (1);
    udn_2412[indx].record3.push_back (l1);
    udn_2412[indx].record3.push_back (l2);
    return label;              
    
  }
  
  unsigned Mesh_modifier::add_quad_to_2412 (unsigned l1, unsigned l2, unsigned l3, unsigned l4) {
    unsigned uci = unv_container.size() - 1;
    auto & udn_2412 = unv_container[uci].udn_2412;
    udn_2412.push_back(Universal_dataset_number_2412());
    auto indx = udn_2412.size() - 1;
    auto last_label = udn_2412[indx-1].record1[0];
    auto label = last_label + 1;    
    udn_2412[indx].record1.push_back (label);
    udn_2412[indx].record1.push_back (44);
    udn_2412[indx].record1.push_back (2);
    udn_2412[indx].record1.push_back (1);
    udn_2412[indx].record1.push_back (7);
    udn_2412[indx].record1.push_back (4);
    udn_2412[indx].record2.push_back (l1);
    udn_2412[indx].record2.push_back (l2);
    udn_2412[indx].record2.push_back (l3);
    udn_2412[indx].record2.push_back (l4);
    
    return label;              

  }  
  
  unsigned Mesh_modifier::add_hexahedron_to_2412 (unsigned l1, unsigned l2, unsigned l3, unsigned l4
                                               ,unsigned l5, unsigned l6, unsigned l7, unsigned l8) {
    unsigned uci = unv_container.size() - 1;
    auto & udn_2412 = unv_container[uci].udn_2412;
    udn_2412.push_back(Universal_dataset_number_2412());
    auto indx = udn_2412.size() - 1;
    auto last_label = udn_2412[indx-1].record1[0];
    auto label = last_label + 1;    
    udn_2412[indx].record1.push_back (label);
    udn_2412[indx].record1.push_back (115);
    udn_2412[indx].record1.push_back (2);
    udn_2412[indx].record1.push_back (1);
    udn_2412[indx].record1.push_back (7);
    udn_2412[indx].record1.push_back (8);
    udn_2412[indx].record2.push_back (l1);
    udn_2412[indx].record2.push_back (l2);
    udn_2412[indx].record2.push_back (l3);
    udn_2412[indx].record2.push_back (l4);
    udn_2412[indx].record2.push_back (l5);
    udn_2412[indx].record2.push_back (l6);
    udn_2412[indx].record2.push_back (l7);
    udn_2412[indx].record2.push_back (l8);    

    return label;              
  }  
  
  

  


}


