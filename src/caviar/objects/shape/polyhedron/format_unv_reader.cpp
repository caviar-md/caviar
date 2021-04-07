
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

#include "caviar/objects/shape/polyhedron/format_unv_reader.h"
#include "caviar/objects/shape/polyhedron/polyhedron.h"
#include "caviar/objects/shape/polyhedron/preprocess.h"
#include "caviar/utility/interpreter_io_headers.h"

#include <iomanip>

namespace caviar {
namespace objects {
namespace shape {
namespace polyhedron {

Format_unv_reader::Format_unv_reader (CAVIAR *fptr) : Pointers{fptr} {}

Format_unv_reader::~Format_unv_reader() {}

 
void Format_unv_reader::read_polyhedron (shape::polyhedron::Polyhedron &p_object, const std::string & filename) {

    std::ifstream ifs;
    ifs.open (filename.c_str());

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
        if (udn_code == 2411) import_udn_2411 (p_object, ifs);
        else if (udn_code == 2412) import_udn_2412 (p_object, ifs);
        else if (udn_code == 2467) import_udn_2467 (p_object, ifs);
        else import_udn_ignore (ifs, udn_code);
        
        in_section = false;
      }
    }
    ifs.close ();
}

void Format_unv_reader::import_udn_ignore (std::ifstream & ifs, int udn_code) {
    std::cout << "Warning: Unsupported udn " << udn_code 
              << ". Ignoring the section.\n";  
    while (true) {
     std::string c;
     ifs >> c;
     if (c == "-1") return;
    }  
}  

  
void Format_unv_reader::import_udn_2411 (shape::polyhedron::Polyhedron &p_object, std::ifstream & ifs) {

  auto & vertex = p_object.vertex;


    while (true) {
      int tmp;
      ifs >> tmp;
      if (tmp == -1) {
        polyhedron::Preprocess p_pre (fptr);
        p_pre.merge_vertices(p_object); 

        return;
      }
      
      int dummy[3];
      ifs >> dummy[0] >> dummy[1] >> dummy[2];
      //std::cout<< tmp << " " <<  dummy[0] << " " << dummy[1] << " " << dummy[2] << "\n";      
     
      double pos[3];      
      ifs >> pos[0] >> pos[1] >> pos[2];
      //std::cout << pos[0] << " " << pos[1] << " " << pos[2] << "\n";       

      node_label.push_back(tmp);

      vertex.push_back (Vector<Real_t> {pos[0], pos[1], pos[2]});

    }
     
}

void Format_unv_reader::import_udn_2412 (shape::polyhedron::Polyhedron &p_object, std::ifstream & ifs) {
  std::cout << "import_udn_2412" << std::endl;

  auto & edges = p_object.edges;
  auto & face = p_object.face;
  auto & face_id = p_object.face_id;
  auto & vertex_map = p_object.vertex_map;

  if (face.size() != face_id.size()) {
    error->all(FC_FILE_LINE_FUNC,"face.size() != face_id.size()");
  }


    while (true) {
      int tmp;
      ifs >> tmp;
      if (tmp == -1) return;


      
      int dummy[5];
      ifs >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> dummy[4];
      
      //std::cout << tmp << " " <<  dummy[0] << " " << dummy[1] << " " << dummy[2] 
      //          << " " <<  dummy[3] << " " << dummy[4] << "\n";




      int FE_Id = dummy[0];
      bool beam_type = (FE_Id==11 || FE_Id==21 || FE_Id==22 || FE_Id==23 || FE_Id==24);
      unsigned num_of_elements = dummy[4];
                            
      if (beam_type) { //beam elements
        int dummy_[2];
        ifs >> dummy_[0] >> dummy_[1] >> dummy_[2];
        //std::cout <<  dummy_[0] << " " << dummy_[1] << " " << dummy_[2]  << "\n";

        
        //int field[num_of_elements];
        std::vector<int> field(num_of_elements);
        for (unsigned int i = 0; i < num_of_elements; ++i)
          ifs >> field[i];

        //for (unsigned int i = 0; i < num_of_elements; ++i)
          //std::cout << field[i] << " ";
          //std::cout << "\n";
         
      } else {
        //int field[num_of_elements];
        std::vector<int> field(num_of_elements);
        for (unsigned int i = 0; i < num_of_elements; ++i) {
          ifs >> field[i];
        }

        // if it is a triangle
        if (num_of_elements==3 ) {

          face_label.push_back(tmp);

          unsigned num1 = field[0];
          unsigned num2 = field[1];
          unsigned num3 = field[2];

          num1 = (vertex_map[num1-1].size()==0)?num1-1 :vertex_map[num1-1][1];  // uses vertex_map instead of vertex
          num2 = (vertex_map[num2-1].size()==0)?num2-1 :vertex_map[num2-1][1];  //
          num3 = (vertex_map[num3-1].size()==0)?num3-1 :vertex_map[num3-1][1];  //


          std::vector<unsigned int> gons{num1, num2, num3};
          face.push_back (gons);      
          face_id.push_back(-1); // invalid face_id;
          unsigned this_face = face.size() - 1;

          std::map<std::vector<unsigned int>,std::vector<unsigned int>>::iterator it_edges; 
#define   ADD_EDGE(NUM1,NUM2)                                    \
          {                                                      \
            std::vector<unsigned int> check_face = {this_face};  \
            std::vector<unsigned int> temp_edge;                 \
            if (NUM1<NUM2) temp_edge = {NUM1,NUM2};              \
            else temp_edge = {NUM2,NUM1};                        \
            it_edges = edges.find(temp_edge);                    \
            if (it_edges == edges.end()) {                       \
              edges.insert (make_pair(temp_edge,check_face));    \
            }  else {                                            \
              it_edges->second.push_back(this_face);             \
            }                                                    \
          }
          ADD_EDGE(num1,num2);
          ADD_EDGE(num2,num3);
          ADD_EDGE(num3,num1);
#undef ADD_EDGE

        // if it is a triangle
        } else if (num_of_elements==4 && FE_Id != 111) {
          error->all(FC_FILE_LINE_FUNC,"(num_of_elements==4 ), not implemented yet");
        }
  
      }

    }


     
}
  
  
void Format_unv_reader::import_udn_2467 (shape::polyhedron::Polyhedron &p_object, std::ifstream & ifs) {

  auto & face_id = p_object.face_id;

  make_label_to_index(face_label, face_label_to_index);


    while (true) {

      int tmp;
      ifs >> tmp;
      if (tmp == -1) return;

      
      int dummy[7];

      ifs >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> dummy[4]
          >> dummy[5] >> dummy[6] ;


        
      int group_name;
      ifs >> group_name;
  
      //std::cout << "group_name: " << group_name << "\n";


      for (int i = 0; i < dummy[6]; ++i) {
        int dummy[4];
        ifs >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] ;        


        int inx = face_label_to_index[dummy[1]];

        //std::cout << dummy[0] << " " << dummy[1] << " " << dummy[2] << " " 
        //          << dummy[3] ;
         


      
        face_id[inx] = group_name; 


      }


    }
}


//============================================
//============================================
//============================================
//============================================


void Format_unv_reader::write_unv (shape::polyhedron::Polyhedron &p_object, const std::string filename) {

    
    std::ofstream ofs;
    ofs.open (filename.c_str());
        
   
    export_udn_2411 (p_object, ofs);
    export_udn_2412 (p_object, ofs);
    export_udn_2467 (p_object, ofs);
    
    ofs.close ();    
}
  
  
    
void Format_unv_reader::export_udn_2411 (shape::polyhedron::Polyhedron &p_object, std::ofstream & ofs) {

  const auto & vertex = p_object.vertex;
  //auto & vertex_map = p_object.vertex_map;

  ofs << std::right << std::setw(6) << "-1" << "\n";               
  ofs << std::right << std::setw(6) << "2411" << "\n";                   
  for (unsigned i = 0; i < vertex.size(); ++i) {
      ofs << std::right 
          << std::setw(10) << i + 1 // Plus one is due to starting label from 1 in UNV file.
          << std::setw(10) << 1
          << std::setw(10) << 1
          << std::setw(10) << 11 << "\n";      
      

      ofs << std::scientific << std::setprecision(16) << std::uppercase
          << std::setw(25) << vertex[i].x 
          << std::setw(25) << vertex[i].y
          << std::setw(25) << vertex[i].z << "\n"; 
      ofs << std::fixed;
           
  }
  ofs << std::right << std::setw(6) << "-1" << "\n";          
    
}
  

void Format_unv_reader::export_udn_2412 (shape::polyhedron::Polyhedron &p_object,std::ofstream & ofs) {

    const auto & face = p_object.face;
    
  
    ofs << std::right << std::setw(6) << "-1" << "\n";               
    ofs << std::right << std::setw(6) << "2412" << "\n";
                       
    for (unsigned i = 0; i < face.size(); ++i) {
      auto num_of_elements = face[i].size();
      int fe_descriptor_id = 0;
      if (num_of_elements==3)
        fe_descriptor_id = 41;
      else if (num_of_elements==4)
        fe_descriptor_id = 44;
      else continue;

      ofs << std::right 
          << std::setw(10) << i + 1 // Plus one is due to starting label from 1 in UNV file.
          << std::setw(10) << fe_descriptor_id
          << std::setw(10) << 2
          << std::setw(10) << 1
          << std::setw(10) << 7
          << std::setw(10) << num_of_elements << "\n";  
              

      
      for (unsigned int j = 0; j < num_of_elements; ++j)     
        ofs << std::setw(10) << face[i][j]+1;      // Plus one is due to starting label from 1 in UNV file.
      ofs << "\n";                    
        
      
    }
    ofs << std::right << std::setw(6) << "-1" << "\n";  
}                  
   
                

void Format_unv_reader::export_udn_2467 (shape::polyhedron::Polyhedron &p_object, std::ofstream & ofs) {

  auto & face_id = p_object.face_id;


  int id_max = -9999;
  for (unsigned i = 0; i < face_id.size(); ++i) {
    if (id_max < face_id[i]) id_max = face_id[i];
  }

  if (id_max < 0) return;

  std::vector<std::vector<int>> id_list(id_max+1);
  for (unsigned i = 0; i < face_id.size(); ++i) {
    if (face_id[i] > -1)
      id_list[face_id[i]].push_back(i);
  }

  ofs << std::right << std::setw(6) << "-1" << "\n";               
  ofs << std::right << std::setw(6) << "2467" << "\n";
  int counter = 0;
  for (unsigned i = 0; i < id_list.size(); ++i) {  

      unsigned num_of_elements = id_list[i].size();
      if (num_of_elements == 0) continue;

      counter++;
      ofs << std::right 
          << std::setw(10) << counter   // just a counter
          << std::setw(10) << 0
          << std::setw(10) << 0
          << std::setw(10) << 0
          << std::setw(10) << 0
          << std::setw(10) << 0
          << std::setw(10) << 0
          << std::setw(10) << num_of_elements << "\n";
                                                                                                
     //group name
      ofs << i << "\n";
      
      for (unsigned j = 0; j < num_of_elements; ++j) {
        ofs << std::setw(10) << 8 
            << std::setw(10) << id_list[i][j] + 1  // Plus one is due to starting label from 1 in UNV file.
            << std::setw(10) << 0
            << std::setw(10) << 0;
        if (j%2==1)
          ofs << "\n";
      }
    
      if (num_of_elements%2 == 1)
        ofs << "\n";
      
    }   
    
    ofs << std::right << std::setw(6) << "-1" << "\n";  
} 


//============================================
//============================================
//============================================
//============================================


void Format_unv_reader::make_label_to_index (const std::vector<int> & u
                                            ,std::vector<int> & v) {

  int max_label = 0;
  for (unsigned i = 0; i < u.size(); ++i) 
    if (max_label < u[i]) 
      max_label = u[i];


  v.resize(max_label + 1, -1);
  for (unsigned i = 0; i < u.size(); ++i) 
    v[u[i]] = i;            
      
}

} //polyhedron
} //shape
} //objects
} // namespace caviar

