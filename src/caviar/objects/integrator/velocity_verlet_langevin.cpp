
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

#include "caviar/objects/integrator/velocity_verlet_langevin.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/force_field.h"


namespace caviar {
namespace objects {
namespace integrator {

Velocity_verlet_langevin::Velocity_verlet_langevin (CAVIAR *fptr) : Integrator{fptr} {
  FC_OBJECT_INITIALIZE_INFO
  integrator_type = 3;
  rnd_generator_x.seed(1);
  rnd_generator_y.seed(2);
  rnd_generator_z.seed(3);
  kb = -1; temperature = -1; kbt = -1; friction = 0.0; dt = -1.0;
  initialized = false;
}

Velocity_verlet_langevin::~Velocity_verlet_langevin(){}

bool Velocity_verlet_langevin::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"temperature")) {
      GET_OR_CHOOSE_A_REAL(temperature,"","")
      if (temperature < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "Temperature have to non-negative."); 
    } else if (string_cmp(t,"friction")) {
      GET_OR_CHOOSE_A_REAL(friction,"","")
      if (friction < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "friction have to non-negative."); 
    } else if (string_cmp(t,"kb")) {
      GET_OR_CHOOSE_A_REAL(kb,"","")
      if (kb < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "kb have to non-negative."); 
    } else if (string_cmp(t,"kbt")) {
      GET_OR_CHOOSE_A_REAL(kbt,"","")
      if (kb < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "kb have to non-negative."); 
    } else if (string_cmp(t,"print_langevin_parameters")) {
      print_langevin_parameters();
    } else if (string_cmp(t,"dt")) {
      GET_OR_CHOOSE_A_REAL(dt,"","")
    } else FC_ERR_UNDEFINED_VAR(t)
  }

  return in_file;
}


void Velocity_verlet_langevin::verify_settings (){
  FC_NULLPTR_CHECK(atom_data)
}


void Velocity_verlet_langevin::print_langevin_parameters() {
  auto &vel = atom_data -> owned.velocity;
  Vector<double> sum_v {0.0,0.0,0.0};
  for (unsigned int i=0; i<vel.size(); i++) {
  
    std:: cout << i << " : " << vel [i] << " , " << vel[i]/dt << std::endl;
    sum_v += vel[i];
  }

  std::cout << "print_langevin_parameters:\n"
            << "a: " << a << " b: " << b << " c: " << c << "\n"
            << "b/(0.5*dt): " << b/(0.5*dt) << "(acceleration dimension)\n" 
            << "c*a: " << c*a << "(position dimension)\n" 
            << "c*b: " << c*b << "(position dimension)\n" << std::flush;
}

void Velocity_verlet_langevin::initialize () {
  if (dt < 0.0)
    error->all(FC_FILE_LINE_FUNC,"expected 'dt' input");    

  a = (2.0 - friction * dt) / (2.0 + friction * dt);
  if (kbt > 0)
    b = std::sqrt(kbt * friction * 0.5 * dt);
  else if( kb >0 && temperature > 0)
    b = std::sqrt(kb  * temperature * friction * 0.5 * dt);  
  else
    error->all(FC_FILE_LINE_FUNC,"expected ('temperature'  & 'kb') or 'kbt' input");
  c = 2.0 * dt / (2.0 + friction * dt);
  initialized = true;
}


// It is introduced as 'Classical Velocity Verlet algorithm' and 
// the Langevin method for it according to
// M. Kröger, 'Langevin dynamics modified Velocity Verlet algorithm'
// in which cited, M. Kröger,
// Models for polymeric and anisotropic liquids
// (Springer, Berlin, 2005)
// However, the 'Classical Velocity Verlet algorithm' is called as another 
// implementaion of Leapfrog at Rapaport 'Art of Molecular Dynamics'. 
//
void Velocity_verlet_langevin::step_part_I () {

  FC_OBJECT_VERIFY_SETTINGS

  if (!initialized) {
    initialize();
  }

  //auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;

  const auto psize = vel.size();
  
  eta_x.resize(psize);
  eta_y.resize(psize);
  eta_z.resize(psize);

#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 


    eta_x[i] = rnd_ndist_x (rnd_generator_x);      
    eta_y[i] = rnd_ndist_y (rnd_generator_y); 
    eta_z[i] = rnd_ndist_z (rnd_generator_z);
    
    const auto eta = Vector<double>{eta_x[i], eta_y[i], eta_z[i]};

    vel [i] += 0.5 * acc [i] * dt + b * eta;


  }
}

void Velocity_verlet_langevin::step_part_II () {
  auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;
  const auto psize = pos.size();

#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 
    pos [i] += vel [i] * c;
  }
}

void Velocity_verlet_langevin::step_part_III () {

  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;

#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<vel.size(); i++) {
    const auto eta = Vector<double>{eta_x[i], eta_y[i], eta_z[i]};
    vel [i] = a * vel [i] + b * eta + 0.5 * acc [i] * dt;

  }
}



} //integrator
} //objects

} // namespace caviar


