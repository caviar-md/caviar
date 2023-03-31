
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

#include "caviar/objects/force_field/granular.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"
#include "caviar/utility/vector.h"
#include "caviar/utility/interpreter_io_headers.h"
#include <cmath>

CAVIAR_NAMESPACE_OPEN

namespace force_field {

Granular::Granular (CAVIAR *fptr) : Force_field{fptr} {
  FC_OBJECT_INITIALIZE_INFO
}

bool Granular::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"cutoff")) {
      GET_OR_CHOOSE_A_REAL(cutoff,"","")
      if (cutoff < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "Force field cutoff have to non-negative.");      
    } else if (string_cmp(t,"force_type")) {
      GET_OR_CHOOSE_A_INT(force_type,"","")
      if (force_type > 1 || force_type < 0)
        error->all (FC_FILE_LINE_FUNC_PARSE, "invalid force_type");      
    } else if (string_cmp(t,"gravity")) {
      GET_OR_CHOOSE_A_REAL_3D_VECTOR(gravity,"","")
    } else if (string_cmp(t,"radius")) {
      GET_A_STDVECTOR_REAL_ELEMENT(radius)
      if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "radius have to be non-negative.");      
    } else if (string_cmp(t,"elastic_coef")) {
      GET_A_STDVECTOR_REAL_ELEMENT(elastic_coef)
      if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Elastic coef. have to be non-negative.");      
    }  else if (string_cmp(t,"dissip_coef")) {
      GET_A_STDVECTOR_REAL_ELEMENT(dissip_coef)
      //if (vector_value < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "Dissipation coef. have to be non-negative.");            
    } else if (string_cmp(t,"set_neighborlist") || string_cmp(t,"neighborlist")) {
      FIND_OBJECT_BY_NAME(neighborlist,it)
      neighborlist = object_container->neighborlist[it->second.index];
    } else if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"set_domain") || string_cmp(t,"domain")) {
      FIND_OBJECT_BY_NAME(domain,it)
      domain = object_container->domain[it->second.index];
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  
  return in_file;
}



void Granular::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)
  FC_NULLPTR_CHECK(neighborlist)
  FC_NULLPTR_CHECK(domain)
}


void Granular::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS

  //const auto &g_pos = atom_data -> ghost.position;
  //const auto &g_vel = atom_data -> ghost.velocity;
  //const auto &g_type = atom_data -> ghost.type;
  //const auto &g_id = atom_data -> ghost.id;

  auto a_radius = atom_data->owned.radius;

  auto cutoff_sq = cutoff * cutoff;
  const auto &nlist = neighborlist -> neighlist;
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned i=0; i<nlist.size (); ++i) {
    const auto &pos_i = atom_data -> owned.position [i];
    const auto &vel_i = atom_data -> owned.velocity [i];
    const auto type_i = atom_data -> owned.type [i];
    const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];
    const auto ela_i = elastic_coef [type_i];
    const auto dis_i = dissip_coef [type_i];
    const auto rad_i = a_radius [type_i];

    if (mass_inv_i > 0)
         atom_data -> owned.acceleration [i] += gravity;
    for (auto j : nlist[i]) {
      //std::cout << "s: " << atom_data -> owned.acceleration .size() << "\n";
      //std::cout << "i: " << i << "\n";
      //std::cout << "j: " << j << std::endl;
      //std::cout << nlist[i].size()  << " " ;    
      bool is_ghost = j >= nlist.size();
      Vector<Real_t> pos_j, vel_j;
      Real_t type_j, mass_inv_j;
      if (is_ghost) {
        j -= nlist.size ();
        //std::cout << "jg: " << j << std::endl;
        pos_j = atom_data->ghost.position [j];
        vel_j = atom_data->ghost.velocity [j];
        type_j = atom_data->ghost.type [j];
      } else {
        pos_j = atom_data->owned.position [j];
        vel_j = atom_data->owned.velocity [j];
        type_j = atom_data->owned.type [j];
      }
      mass_inv_j = atom_data->owned.mass_inv [ type_j ];

      const auto dr = pos_j - pos_i;

      auto r_sq = dr*dr;

      if (r_sq > cutoff_sq) continue;
      const auto ela_j = elastic_coef [type_j];
      const auto dis_j = dissip_coef [type_j];
      const auto rad_j = a_radius [type_j];

      auto rr = std::sqrt(r_sq);
      auto xi = rad_i+rad_j-rr;

      if (xi > 0.0) {

        auto Y = ela_i*ela_j/(ela_i+ela_j);
        auto A = -0.5*(dis_i+dis_j);
        auto reff = (rad_i*rad_j)/(rad_i+rad_j);
        auto dv = vel_i - vel_j;

        auto rr_rez = 1.0/rr;
        auto e = rr_rez * dr;
        auto xidot = - e * dv;

        double fn;
        switch (force_type){
          case 0: default:
          fn =Y*std::sqrt(reff)*(+xi+A*xidot);                    //linear-dashpot        
          break;

          case 1:
          fn = std::sqrt(xi)*Y*std::sqrt(reff)*(+xi+A*xidot);  //visco-elastic
          break;
        }
        if (fn < 0.0) fn = 0.0;
        auto force = -fn * e;



        if (mass_inv_i > 0)
          atom_data -> owned.acceleration [i] += force * mass_inv_i;

        if (!is_ghost) {
          if (mass_inv_j > 0) {
#ifdef CAVIAR_WITH_OPENMP        
#pragma omp atomic 
          atom_data -> owned.acceleration [j].x -= force.x * mass_inv_j;   
#pragma omp atomic
          atom_data -> owned.acceleration [j].y -= force.y * mass_inv_j;   
#pragma omp atomic 
          atom_data -> owned.acceleration [j].z -= force.z * mass_inv_j;  
#else
          atom_data -> owned.acceleration [j] -= force * mass_inv_j;
#endif    
            
           // atom_data -> owned.acceleration [j] -= force * mass_inv_j;
          }
        }

      }       

    }   

  }
}

} //force_field

CAVIAR_NAMESPACE_CLOSE

