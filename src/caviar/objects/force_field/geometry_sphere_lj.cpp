
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

#include "caviar/objects/force_field/geometry_sphere_lj.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"
#include <string>
#include <cmath>
#include <fstream>

namespace caviar {
namespace objects {
namespace force_field {

Geometry_sphere_lj::Geometry_sphere_lj (CAVIAR *fptr) : Force_field {fptr}
{ 
  FC_OBJECT_INITIALIZE_INFO
  wca = false;
  cutoff_list_activated = false;
  force_coef = 1.0;
  center = caviar::Vector<double> (0,0,0);  
  inside = true;
  radius = 1.0;
}

Geometry_sphere_lj::~Geometry_sphere_lj () { 

}

bool Geometry_sphere_lj::read (class caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"cutoff")) {
      GET_OR_CHOOSE_A_REAL(cutoff,"","")
      if (cutoff < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "Force field cutoff have to non-negative.");    
    } else if (string_cmp(t,"cutoff_list")) {
      GET_A_STDVECTOR_REAL_ELEMENT(cutoff_list)
      cutoff_list_activated = true;
      if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Epsilon have to be non-negative.");      
    } else if (string_cmp(t,"wca")) {
      wca = true;
      cutoff_list_activated = true;
    } else if (string_cmp(t,"center")) {
      GET_OR_CHOOSE_A_REAL_3D_VECTOR(center,"","")
    } else if (string_cmp(t,"radius")) {
      GET_OR_CHOOSE_A_REAL(radius,"","")
      if (radius < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "radius have to be non-negative.");      
    } else if (string_cmp(t,"outside")) {
      inside = false;
    } else if (string_cmp(t,"inside")) {
      inside = true;
    } else if (string_cmp(t,"force_coef")) {
      GET_OR_CHOOSE_A_REAL(force_coef,"","")
    } else if (string_cmp(t,"epsilon_atom")) {
      GET_A_STDVECTOR_REAL_ELEMENT(epsilon_atom)
    } else if (string_cmp(t,"epsilon_wall")) {
      GET_OR_CHOOSE_A_REAL(epsilon_wall,"","")
    } else if (string_cmp(t,"sigma_atom")) {
      GET_A_STDVECTOR_REAL_ELEMENT(sigma_atom)
    } else if (string_cmp(t,"sigma_wall")) {
      GET_OR_CHOOSE_A_REAL(sigma_wall,"","")
    } else if (string_cmp(t,"set_neighborlist") || string_cmp(t,"neighborlist")) {
      FIND_OBJECT_BY_NAME(neighborlist,it)
      neighborlist = object_container->neighborlist[it->second.index];
    } else if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  
  return in_file;
}


void Geometry_sphere_lj::verify_settings() {
    FC_NULLPTR_CHECK(atom_data)

  
  auto epsilon_atom_size = epsilon_atom.size();
  auto sigma_atom_size   = sigma_atom.size();
  auto atom_max_size = (epsilon_atom_size > sigma_atom_size ? epsilon_atom_size : sigma_atom_size);
  if (epsilon_atom_size != sigma_atom_size) {
    output->warning("Geometry_lj:: (epsilon_atom_size != sigma_atom_size)");
    if (epsilon_atom_size != atom_max_size) {
      epsilon_atom.resize(atom_max_size,0);
    } else {
      sigma_atom.resize(atom_max_size,0);
    }
  }


  unsigned atom_data_type_max = 0;
  for (auto && t : atom_data->owned.type) {
    if (atom_data_type_max < t) atom_data_type_max = t;
  }


  if (atom_max_size < atom_data_type_max + 1) {
    output->warning("Geometry_lj:: (atom_max_size < atom_data_type_max + 1)");
    epsilon_atom.resize(atom_data_type_max + 1,0);
      sigma_atom.resize(atom_data_type_max + 1,0);
    atom_max_size = atom_data_type_max + 1;
  }

  epsilon.resize(atom_max_size,0);
  sigma.resize(atom_max_size,0);


  
    
    for (unsigned int j = 0; j < atom_max_size; ++j) {
      auto sigma_ij   = 0.5*(sigma_wall + sigma_atom[j]);
      auto epsilon_ij = std::sqrt(epsilon_wall * epsilon_atom[j]);
      sigma[j] = sigma_ij;
      epsilon[j] = epsilon_ij;
    }
  

  //Week-Chandler-Anderson (WCA) potential activated.
  if (wca)  {
    auto cut_coef   = std::pow(2.0, 1.0/6.0);   // ONLY for LJ 6-12

    cutoff_list.resize(atom_max_size,0); 
    
    for (unsigned int j = 0;   j < atom_max_size; ++j) {
      auto cut = cut_coef * sigma[j];
      std::cout << cut << "\n";
      cutoff_list[j] = cut;
    }    
  } 
}


void Geometry_sphere_lj::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS


  const auto &pos = atom_data -> owned.position;
  auto &acc = atom_data -> owned.acceleration;
  
  
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
    for (unsigned int i=0;i<pos.size();++i) {
      const auto type_i = atom_data -> owned.type [i] ;
      const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];   
                 
            
      double c = cutoff;    
     
      if (cutoff_list_activated) c = cutoff_list[type_i];
      
      auto cp = (pos[i]-center);
      auto cp_sq = cp*cp;
      auto cp_abs = std::sqrt(cp_sq);
      auto cp_norm = cp / cp_abs;
      auto force_coef = -1.0;
      
      //double compression;
      
      if (inside)
      {
          force_coef = 1.0;
      }
              
      
      auto d = radius - cp_abs; // distance between wall's atom and atom.
      
      bool is_in_contact = (radius - cp_abs < c);
      
      // if distance to the wall is less than cutoff.
      if (is_in_contact) { 


     

        auto dr = cp_norm * d;
        auto dr_sq = dr*dr;

        auto eps_ij   = epsilon [type_i];
        auto sigma_ij =  sigma [type_i];

        auto r_c_sq_inv    = 1/(c*c);
        auto rho_c_sq_inv  = sigma_ij*sigma_ij*r_c_sq_inv;
        auto rho_c_6_inv   = rho_c_sq_inv*rho_c_sq_inv*rho_c_sq_inv;
        auto rho_c_12_inv  = rho_c_6_inv*rho_c_6_inv;
      
        auto dr_sq_inv     = 1/dr_sq;
        auto rho_sq_inv    = sigma_ij*sigma_ij*dr_sq_inv;
        auto rho_6_inv     = rho_sq_inv*rho_sq_inv*rho_sq_inv;
        auto rho_12_inv    = rho_6_inv*rho_6_inv;

        auto force = force_coef 
                   * 4*eps_ij*(-12*rho_12_inv*dr_sq_inv + 6*rho_6_inv*dr_sq_inv + 
                   +12*rho_c_12_inv*r_c_sq_inv - 6*rho_c_6_inv*r_c_sq_inv   ) * dr;
        //std::cout << "f : " << force << "\n";
        

        acc[i] += force * mass_inv_i;

      }
      
    }

    
    
}

} //force_field
} //objects
} // namespace caviar

