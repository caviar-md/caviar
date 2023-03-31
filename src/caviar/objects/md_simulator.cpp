
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

#include "caviar/objects/md_simulator.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/force_field.h"
#include "caviar/objects/constraint.h"
#include "caviar/objects/writer.h"
#include "caviar/objects/unique/time_function.h"
#include "caviar/objects/unique/time_function_3d.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/utility/time_utility.h"
#include "caviar/interpreter/communicator.h"

namespace caviar {



Md_simulator::Md_simulator (CAVIAR *fptr) : Pointers{fptr},
    atom_data{nullptr},
    initialized{false} {
      
  current_step = 0;
  use_time = false;
  use_step = false;
  
  dt = -1.0;
  
  auto &L = langevin_param;
  L.rnd_generator_x.seed(1);
  L.rnd_generator_y.seed(2);
  L.rnd_generator_z.seed(3);
  L.kb = -1; L.temperature = -1; L.kbt = -1; L.friction = 0.0; 
  
  FC_OBJECT_INITIALIZE
}

Md_simulator::~Md_simulator () {}


bool Md_simulator::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"dt")) {
      GET_OR_CHOOSE_A_REAL(dt,"","")
    } else if (string_cmp(t,"initial_time")) {
      GET_OR_CHOOSE_A_REAL(initial_time,"","")
      use_time = true;
    } else if (string_cmp(t,"final_time")) {
      GET_OR_CHOOSE_A_REAL(final_time,"","")
      use_time = true;
    } else if (string_cmp(t,"initial_step")) {
      GET_OR_CHOOSE_A_INT(initial_step,"","")
      use_step = true;
    } else if (string_cmp(t,"final_step")) {
      GET_OR_CHOOSE_A_INT(final_step,"","")
      use_step = true;
    } else if (string_cmp(t,"run")) {
      run();
    } else if (string_cmp(t,"add_constraint") || string_cmp(t,"constraint")) {
      FIND_OBJECT_BY_NAME(constraint,it)
      constraint.push_back(object_container->constraint[it->second.index]);
    } else if (string_cmp(t,"remove_constraint")) {
#define FC_REMOVE_OBJECT(name) \
      FIND_OBJECT_BY_NAME(name,it) \
      bool found = false; \
      for (unsigned int i = 0 ; i < name .size(); ++i ) \
      { \
          if (name [i] == object_container-> name [it->second.index]) \
          { \
              name .erase(name .begin()+i); \
              found = true; \
              break; \
          } \
      } \
      if (!found) error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown " #name " to remove.");
      FC_REMOVE_OBJECT(constraint) 
    } else if (string_cmp(t,"add_writer") || string_cmp(t,"writer")) {
      FIND_OBJECT_BY_NAME(writer,it)
      writer.push_back(object_container->writer[it->second.index]);
    } else if (string_cmp(t,"add_neighborlist") || string_cmp(t,"neighborlist")) {
      FIND_OBJECT_BY_NAME(neighborlist,it)
      neighborlist.push_back(object_container->neighborlist[it->second.index]);
    } else if (string_cmp(t,"add_force_field") || string_cmp(t,"force_field")) {
      FIND_OBJECT_BY_NAME(force_field,it)
      force_field.push_back(object_container->force_field[it->second.index]);
    } else if (string_cmp(t,"remove_writer")) {
      FC_REMOVE_OBJECT(writer) 
    } else if (string_cmp(t,"remove_neighborlist")) {
      FC_REMOVE_OBJECT(neighborlist)       
    } else if (string_cmp(t,"remove_force_field")) {
      FC_REMOVE_OBJECT(force_field) 
    } else  if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"add_time_function") || string_cmp(t,"time_function")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,time_function)
      unique::Time_function *a = dynamic_cast<unique::Time_function *>(object_container->unique[it->second.index]);
      time_function.push_back(a);            
    } else if (string_cmp(t,"add_time_function_3d") || string_cmp(t,"time_function_3d")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,time_function_3d)
      unique::Time_function_3d *a = dynamic_cast<unique::Time_function_3d *>(object_container->unique[it->second.index]);
      time_function_3d.push_back(a);            
    } else  if (string_cmp(t,"step")) {
      //int i=0;
      //GET_OR_CHOOSE_A_INT(i,"","")
      ++current_step;
      step ();
      continue;
    } else if (string_cmp(t,"integrator_type")) {
      
      auto t2 = parser->get_val_token();
      auto ts = t2.string_value;
      if (string_cmp(ts,"leap_frog"))
        integrator_type = Integrator_t::Leap_frog;
      else if (string_cmp(ts,"velocity_verlet"))
        integrator_type = Integrator_t::Velocity_verlet;
      else if (string_cmp(ts,"velocity_verlet_langevin"))
        integrator_type = Integrator_t::Velocity_verlet_langevin;
      else error->all (FC_FILE_LINE_FUNC_PARSE, "invalid or non-implemented integrator");
      
    } else if (string_cmp(t,"temperature")) {
      GET_OR_CHOOSE_A_REAL(langevin_param.temperature,"","")
      if (langevin_param.temperature < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "Temperature have to non-negative."); 
    } else if (string_cmp(t,"friction")) {
      GET_OR_CHOOSE_A_REAL(langevin_param.friction,"","")
      if (langevin_param.friction < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "friction have to non-negative."); 
    } else if (string_cmp(t,"kb")) {
      GET_OR_CHOOSE_A_REAL(langevin_param.kb,"","")
      if (langevin_param.kb < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "kb have to non-negative."); 
    } else if (string_cmp(t,"kbt")) {
      GET_OR_CHOOSE_A_REAL(langevin_param.kbt,"","")
      if (langevin_param.kb < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "kbt have to non-negative."); 
    } else if (string_cmp(t,"dt")) {
      GET_OR_CHOOSE_A_REAL(dt,"","")
    }
    else if (read_base_class_commands(parser)) {
    } else error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
  }
  return in_file;
}

bool Md_simulator::run () {
  output->info ("MD started.");
#ifdef CAVIAR_WITH_OPENMP
  output->info("CAVIAR is started with OpenMP. Max number of threads is " + std::to_string(omp_get_max_threads()));
#endif  
  verify_settings ();
  
  double wall0 = get_wall_time();
  //double cpu0  = get_cpu_time();
      
  //t_start = clock();

  for (current_step = initial_step; current_step < final_step; ++current_step) 
    step();

  double wall1 = get_wall_time();
  //double cpu1  = get_cpu_time();
      
  cleanup(); 

  output->info ("MD finished.");

  //t_end = clock();
  //double simulation_time =  ( (float)t_end - (float)t_start ) / CLOCKS_PER_SEC;


  // #if defined (CAVIAR_WITH_MPI) 
  //   std::string s = "MPI process number" + std::to_string(comm->me) + " : ";
  // #else
  //   std::string s = "";
  // #endif  
  
  std::string s = "md_simulation run call:" ;
  s += "Wall Time = " + std::to_string(wall1 - wall0) + " sec. , ";  
  //s += "CPU Time  = " + std::to_string(cpu1  - cpu0) + " sec. ";                      
   
  output->info (s);
  
  return true; //WARNING
}

void Md_simulator::verify_settings () {
  FC_NULLPTR_CHECK(atom_data)
  if (use_time && use_step)
    error->all ("simulator md: verify_settings: cannot use both 'time' and 'step'.");
  if (!(use_time || use_step))
    error->all ("simulator md: verify_settings: please assign 'time' or 'step' variables.");
}


bool Md_simulator::read_base_class_commands (caviar::interpreter::Parser *) {
  bool command_called = false;
/*
  parser -> keep_current_token();
  auto token = parser -> get_val_token();
  auto t = token.string_value;
  if (string_cmp(t,"output_total_force") ) {
    command_called = true;
    output_total_force();
  } else if (string_cmp(t,"output_total_energy") ) {
    command_called = true;
    output_total_energy();
  }
*/
  return command_called;
}

void Md_simulator::initialize () {
  
  
  FC_NULLPTR_CHECK(atom_data)
  
  if (dt < 0.0)
    error->all(FC_FILE_LINE_FUNC,"expected 'dt' input");    
  
  switch(integrator_type) {

    //case Integrator_t::Verlet : 
    //error->all(FC_FILE_LINE_FUNC,"not implemented");
    //integrate_verlet();
    
    break;
    
    default:
    case Integrator_t::Velocity_verlet: 
    //integrate_velocity_verlet();
    break;
    
    
    case Integrator_t::Leap_frog :   
    //integrate_leap_frog();  
    break;     
    
    
    case Integrator_t::Velocity_verlet_langevin :
      
    auto &L = langevin_param;

    L.a = (2.0 - L.friction * dt) / (2.0 + L.friction * dt);
    if (L.kbt > 0)
      L.b = std::sqrt(L.kbt * L.friction * 0.5 * dt);
    else if( L.kb >0 && L.temperature > 0)
      L.b = std::sqrt(L.kb  * L.temperature * L.friction * 0.5 * dt);  
    else
      error->all(FC_FILE_LINE_FUNC,"expected ('temperature'  & 'kb') or 'kbt' input");
    L.c = 2.0 * dt / (2.0 + L.friction * dt);
  
    // printing the parameters for langevin

    std::cout << "langevin_parameters:\n"
              << "a: " << L.a << " b: " << L.b << " c: " << L.c << "\n"
              << "b/(0.5*dt): " << L.b/(0.5*dt) << "(acceleration dimension)\n" 
              << "c*a: " << L.c*L.a << "(position dimension)\n" 
              << "c*b: " << L.c*L.b << "(position dimension)\n\n" << std::flush;
            
    
    break;
  }
  
  initialized = true;
}

void Md_simulator::re_calculate_acc() {
  // Here there's 'r(t)' stored in 'atom_data->owned.position'
  // First calculate 'a(t)' by using 'r(t)'
  bool make_list = boundary_condition (); 


  if (make_list) 
    for (auto&& n: neighborlist) n->build_neighlist ();
  else {
    for (auto&& n: neighborlist)
      make_list = make_list || n->rebuild_neighlist ();
    if (make_list) 
      for (auto&& n: neighborlist) n->build_neighlist ();
  }

  atom_data -> reset_owned_acceleration();

  for (auto&& f : force_field)
    f -> calculate_acceleration ();

   
}

void Md_simulator::step (int64_t i) {
  current_step = i;
  step();  
}

void Md_simulator::step () {

  FC_OBJECT_VERIFY_SETTINGS

  if (!initialized) {
    initialize();

    for (auto&& c : constraint)
      c -> verify_settings ();

    time = dt*initial_step;

    for (auto&& tf: time_function) tf->update_time_variable (time);    
    for (auto&& tf: time_function_3d) tf->update_time_variable (time); 
    
    setup();
  }

  atom_data->record_owned_old_data();

  
  
  
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  
  /*
  static int my_mpi_rank = comm->me;

  if (my_mpi_rank==0) {

    integrator -> step_part_I ();

    for (auto&& c : constraint)
      c -> step_part_I (i);

    integrator -> step_part_II ();

    for (auto&& c : constraint)
      c -> step_part_II (i);

    bool make_list = boundary_condition (); 

    if (make_list) 
      for (auto&& n: neighborlist) n->build_neighlist ();
    else {
      for (auto&& n: neighborlist)
        make_list = make_list || n->rebuild_neighlist ();
      if (make_list) 
        for (auto&& n: neighborlist) n->build_neighlist ();
    }

  }

  atom_data -> synch_owned_data(0);

  atom_data -> reset_owned_acceleration();

  for (auto&& f : force_field)
    f -> calculate_acceleration ();

  MPI_Barrier(mpi_comm);

  if (my_mpi_rank==0) {

    integrator -> step_part_III ();

    for (auto&& c : constraint)
      c -> step_part_III (i);
  }

  for (auto&& w : writer) 
    w -> write (i, time);
  */
  
#else

  
  switch(integrator_type) {
    
    default:
    case Integrator_t::Velocity_verlet: 
    integrate_velocity_verlet();
    break;
    
    
    case Integrator_t::Leap_frog :   
    integrate_leap_frog();  
    break;     
    
    
    case Integrator_t::Velocity_verlet_langevin :       
    integrate_velocity_verlet_langevin();      
    break;
  }

#endif

  time += dt;
  for (auto&& tf: time_function) tf->update_time_variable (time);
  for (auto&& tf: time_function_3d) tf->update_time_variable (time);
}


bool Md_simulator::boundary_condition () {
  bool result = false;

  // the result only in MPI case can
  // become true, due to exchanging a particle between domains
  result = atom_data -> exchange_owned (); 

  atom_data -> exchange_ghost ();
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#elif defined(CAVIAR_WITH_MPI)
  MPI_Allreduce (MPI::IN_PLACE, &result, 1, MPI::BOOL, MPI::LOR, mpi_comm);
#endif
  return result;
}

void Md_simulator::setup () {
  t_start = clock();
  if (use_time) {
    initial_step = (int) initial_time/dt;
    final_step = (int) final_time/dt;
  }

  FC_NULLPTR_CHECK(atom_data)
  atom_data->record_owned_old_data();

  if (neighborlist.size()==0) output->warning("Md_simulator::setup: neighborlist.size() = 0");
  if (force_field.size()==0) output->warning("Md_simulator::setup: force_field.size() = 0");

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

  atom_data->synch_owned_data(0);

#endif 


  
  for (auto&& n: neighborlist) n->init ();

  atom_data -> exchange_owned (); // isn't neccesary

  atom_data -> exchange_ghost ();

  for (auto&& n: neighborlist) n->build_neighlist ();

  atom_data -> reset_owned_acceleration();

  for (auto&& f : force_field) 
    f -> calculate_acceleration ();

  for (auto&& w : writer) 
    w -> write (0, 0.0);


}



void Md_simulator::cleanup () {

}


void Md_simulator::integrate_velocity_verlet () {
    
  auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;
  auto &acc_old = atom_data -> owned.acceleration_old; //velocity verlet
  //auto &pos_old = atom_data -> owned.position_old; //verlet - also shake and m_shake uses this
  //const auto two_dt_inv = 1.0/(2.0*dt);//verlet
  //-------------------------------
  //-------------------------------
  
  ///XXX
  atom_data->record_owned_acceleration_old = true;
  if (acc_old.size() != acc.size()) acc_old.resize(acc.size());
  //XXX
  
  
  //re_calculate_acc(); // use r(t) to calculate a(t). The forces has to be velocity independent
                        // It has been calculated in previous step
  
  auto psize = pos.size();
  
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 
    
    pos [i] += vel [i] * dt  +  0.5 * acc[i] * dt * dt; // r(t+dt) = r(t) + v(t)*dt + 1/2 * a(t) * dt^2

  }
  
  for (auto&& c : constraint) c -> apply_on_position (current_step); 
  
  re_calculate_acc(); // use r(t+dt) to calculate a(t+dt)
  
  for (auto&& c : constraint) c -> apply_on_acceleration (current_step); 
  
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif    
  for (unsigned int i=0; i<vel.size(); i++) {
    
    vel [i] += 0.5 * dt * (acc [i] + acc_old[i]); // v(t+dt) = v(t) + ( a(t+dt) + a(t) ) * dt / 2
  }
  
  for (auto&& c : constraint) c -> apply_on_velocity (current_step); 
  

  for (auto&& c : constraint) c -> apply (current_step);

  for (auto&& w : writer) w -> write (current_step, time); // pos = r(t+dt) , vel = v(t+dt)

}


void Md_simulator::integrate_velocity_verlet_langevin () {
  auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;
  //auto &acc_old = atom_data -> owned.acceleration_old; //velocity verlet
  //auto &pos_old = atom_data -> owned.position_old; //verlet - also shake and m_shake uses this
  //const auto two_dt_inv = 1.0/(2.0*dt);//verlet
  

  auto &L = langevin_param;

  const auto psize = vel.size();
  
  L.eta_x.resize(psize);
  L.eta_y.resize(psize);
  L.eta_z.resize(psize);

  for (unsigned int i=0; i<psize; i++) { 
    L.eta_x[i] = L.rnd_ndist_x (L.rnd_generator_x);      
    L.eta_y[i] = L.rnd_ndist_y (L.rnd_generator_y); 
    L.eta_z[i] = L.rnd_ndist_z (L.rnd_generator_z);    
  }
  
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) {
    
    const auto eta = Vector<double>{L.eta_x[i], L.eta_y[i], L.eta_z[i]};

    vel [i] += 0.5 * acc [i] * dt + L.b * eta;
    //std::cout << "acc " << acc[i] << ",vel " << vel[i] << "\n";

  }
  
  for (auto&& c : constraint) c -> apply_on_velocity (current_step); 

  #ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 
    //std::cout << "b " << pos[i] ;
    pos [i] += vel [i] * L.c;
    //std::cout << "a " << pos[i] << "\n";
  }
  
  for (auto&& c : constraint) c -> apply_on_position (current_step); 
  
  re_calculate_acc(); // use r(t+dt) to calculate a(t+dt)
  
  for (auto&& c : constraint) c -> apply_on_acceleration (current_step); 
  
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<vel.size(); i++) {
    const auto eta = Vector<double>{L.eta_x[i], L.eta_y[i], L.eta_z[i]};
    vel [i] = L.a * vel [i] + L.b * eta + 0.5 * acc [i] * dt;
    //std::cout << "2 acc " << acc[i] << ",vel " << vel[i] << "\n";

  }
  
  for (auto&& c : constraint) c -> apply (current_step);
  
  for (auto&& w : writer) w -> write (current_step, time); // pos = r(t+dt) , vel = v(t+dt)
  
}


void Md_simulator::integrate_leap_frog () {
  auto &pos = atom_data -> owned.position;
  auto &vel = atom_data -> owned.velocity;
  auto &acc = atom_data -> owned.acceleration;
  //auto &acc_old = atom_data -> owned.acceleration_old; //velocity verlet
  //auto &pos_old = atom_data -> owned.position_old; //verlet - also shake and m_shake uses this
  //const auto two_dt_inv = 1.0/(2.0*dt);//verlet
  
  //-------------------------------
  //-------------------------------
  //re_calculate_acc(); // use r(t) to calculate a(t). The forces has to be velocity independent
                        // It has been calculated in previous step
  auto psize = pos.size();
    
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 
    vel [i] += 0.5 * dt * acc[i]; // v (t + dt/2) = v (t) + (dt/2) a (t) 
    pos [i] += dt * vel[i];       // r (t + dt) = r (t) + dt * v (t + dt/2)      
  }
  
  for (auto&& c : constraint) c -> apply_on_position (current_step); 
  
  re_calculate_acc(); // use r(t+dt) to calculate a(t+dt)
  
  for (auto&& c : constraint) c -> apply_on_acceleration (current_step); 
   
#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<vel.size(); i++) {
    vel [i] += 0.5 * dt * acc [i]; // v (t + dt) = v (t + dt/2) + (dt/2) a (t + dt)
  }     

  for (auto&& c : constraint) c -> apply_on_velocity (current_step);

  for (auto&& c : constraint) c -> apply (current_step);
  
  for (auto&& w : writer) w -> write (current_step, time); // pos = r(t+dt) , vel = v(t+dt)
  
}

/*
void Md_simulator::integrate_verlet () {
{
  error->all(FILE_LINE_FUNC,"not implemented");
  
  if (pos_old.size() != pos.size()) pos_old.resize(pos.size());
  
  atom_data->record_owned_position_old = true;

  re_calculate_acc(); // use r(t) to calculate a(t). The forces has to be velocity independent


#ifdef CAVIAR_WITH_OPENMP
  #pragma omp parallel for
#endif
  for (unsigned int i=0; i<psize; i++) { 
    
    pos [i] = (2.0*pos[i]) - pos_old[i]  + acc [i] * dt * dt;  // r (t + dt) = 2*r(t) − r (t − dt) + a(t) * dt^2 

    
    // XXX note that this is 'v(t)' not 'v(t+dt)'
    vel [i] = domain->fix_distance(pos[i] - pos_old[i]) * two_dt_inv; // v(t) = (r (t + dt) - r (t - dt) )/ (2*dt) 

  }              
  
  // velocity related constraints are useless here
  for (auto&& c : constraint) c -> apply (i); // use unconstrained r(t+dt) and change it to constrained value
  
  
  for (auto&& w : writer) w -> write (i, time); // pos = r(t+dt) , vel = v(t)      
  
}*/

}

} // namespace caviar

