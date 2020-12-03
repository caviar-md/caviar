
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

#include "caviar/objects/atom_data.h"
#include "caviar/interpreter/communicator.h"
#include "caviar/objects/domain.h"
#include "caviar/objects/force_field.h"

#include "caviar/utility/interpreter_io_headers.h"
#include <algorithm>
#include <random>

namespace caviar {

namespace objects {

constexpr auto expected_imbalance_factor = 1.1;

Atom_data::Atom_data (CAVIAR *fptr) : Pointers{fptr}, 
    num_local_atoms{0},
    num_total_atoms{0}, num_atom_types{0},
    synch_owned_data_bcast_details{true},
    ghost_cutoff{0}, domain{nullptr}, cell_list{nullptr} {

  FC_OBJECT_INITIALIZE
  record_owned_position_old = false;
  record_owned_velocity_old = false;
  record_owned_acceleration_old = false;
  make_ghost_velocity = false;
  cutoff_extra = 0.0;
  k_b = -1.0;
  n_r_df = -1;
}

Atom_data::~Atom_data () {
  
}

void Atom_data::reset_owned_acceleration () {
  for (auto &&i : owned.acceleration) {i.x=0.0; i.y=0.0; i.z = 0.0;}
}


void Atom_data::verify_settings () {
  
}

int Atom_data::get_n_r_df() {
  // look at page 114.
  // Philippe H. Hunenberger, Adv. Polym. Sci. (2005) 173:105–149  :
  // N_r = 0 in the presence of stochastic.
  // N_r = 3 under periodic boundary conditions.
  // N_r = 6 under vacuum boundary conditions.
  if (n_r_df == -1) {
    if (stochastic_force_present) n_r_df = 0;
    else {
      auto bc = domain->boundary_condition;
      auto bc_sum = bc.x + bc.y + bc.z;
      if (bc_sum==0) n_r_df = 6;
      else if (bc_sum==3) n_r_df = 3;

      else 
        // I have made this of myself so that it satisfies bc_sum==0,3 conditions.
        // Use it with care.
        n_r_df = 6 - bc_sum; 
    }
  } 

  return n_r_df;
}

int Atom_data::degree_of_freedoms () {
  auto df = 3*owned.position.size();

  auto sum_of_bonds = 0;
  for (auto && i :   owned.atomic_bond_vector)
    sum_of_bonds += i.size();

  auto sum_of_angles = 0;
  for (auto && i : owned.atomic_angle_vector)
    sum_of_angles += i.size();


  return df - sum_of_bonds - sum_of_angles - get_n_r_df();
}


Vector<Real_t> Atom_data::owned_position_cm () {
  Vector<Real_t> p_cm {0.0,0.0,0.0};
  double mass_sum = 0.0;
  auto p_size = owned.position.size();
  for (unsigned int i = 0; i < p_size; ++i) {
    auto type_i = owned.type[i];
    auto mass_i = owned.mass[type_i];
    mass_sum += mass_i;
    p_cm += owned.position[i] * mass_i;
  }
  p_cm = p_cm / mass_sum;
  return p_cm;
}

Vector<Real_t> Atom_data::owned_velocity_cm () {
  Vector<Real_t> v_cm {0.0,0.0,0.0};
  double mass_sum = 0.0;
  auto p_size = owned.velocity.size();
  for (unsigned int i = 0; i < p_size; ++i) {
    auto type_i = owned.type[i];
    auto mass_i = owned.mass[type_i];
    mass_sum += mass_i;
    v_cm += owned.velocity[i] * mass_i;
  }
  v_cm = v_cm / mass_sum;
  return v_cm;
}

Vector<Real_t> Atom_data::owned_angular_momentum_cm () {
  Vector<Real_t> am_cm {0.0,0.0,0.0};
  error->all(FC_FILE_LINE_FUNC,"not implemented.");
    return am_cm;
}

std::vector<std::vector<Real_t> > Atom_data::owned_inertia_tensor_cm () {
  std::vector<std::vector<Real_t> > it_cm (3, std::vector<Real_t> (3, 0.0));
  owned_position_cm();
  error->all(FC_FILE_LINE_FUNC,"not implemented.");
  return it_cm;
}
  

void Atom_data::record_owned_old_data () {
  if (record_owned_position_old) owned.position_old = owned.position;
  if (record_owned_velocity_old) owned.velocity_old = owned.velocity;
  if (record_owned_acceleration_old) owned.acceleration_old = owned.acceleration;
}

unsigned int Atom_data::get_global_id () {
#ifdef CAVIAR_WITH_MPI
  MPI_Barrier (mpi_comm); // does it have to be here??
  MPI_Allreduce (&num_local_atoms, &num_total_atoms, 1, MPI_UNSIGNED, MPI_SUM, mpi_comm);
#else
  num_total_atoms = owned.position.size();
#endif
  return num_total_atoms;
}

void Atom_data::set_num_total_atoms (GlobalID_t n) {
  num_total_atoms = n;
  num_local_atoms_est = n * expected_imbalance_factor / comm->nprocs;
}

void Atom_data::reserve_owned_vectors () {
  owned.id.reserve (num_local_atoms_est);
  owned.charge.reserve (num_local_atoms_est);
  owned.position.reserve (num_local_atoms_est);
  owned.velocity.reserve (num_local_atoms_est);
  owned.acceleration.reserve (num_local_atoms_est);
}

bool Atom_data::position_inside_local_domain(const Vector<double> &pos) {
  if (domain==nullptr) error->all(FC_FILE_LINE_FUNC,"domain = nullptr");

  if (pos.x >= domain->lower_local.x && pos.x < domain->upper_local.x &&
      pos.y >= domain->lower_local.y && pos.y < domain->upper_local.y &&
      pos.z >= domain->lower_local.z && pos.z < domain->upper_local.z) 
    return true;
  
  return false;
}


// XXX note that any new vector addition to this function should be deleted in 
// 'remove_atom()' functions.
bool Atom_data::add_atom (GlobalID_t id,
                          AtomType_t type,
                          const Vector<Real_t> &pos,
                          const Vector<Real_t> &vel,
                          const std::vector<Real_t> &,
                          const std::vector<int> &) {

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  const auto me = domain->me;
  if (me==0) {
#endif
  if (position_inside_local_domain (pos)) {
    owned.type.emplace_back (type);
    owned.position.emplace_back (pos);
    owned.id.emplace_back ( id );
    owned.velocity.emplace_back (vel);
    owned.acceleration.emplace_back (0,0,0);
    owned.msd_domain_cross.emplace_back(0,0,0);
    ++num_local_atoms;
    return true;
  }
  else return false;
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  }
#endif
  return false;
}

bool Atom_data::add_masses (unsigned int type, Real_t m) {
  if (type + 1 > owned.mass.size()) {
    owned.mass.resize (type + 1);
    owned.mass_inv.resize (type + 1);
  }
  owned.mass[type] = m;
  if (m==0.0) owned.mass_inv[type] = 0.0;
  else owned.mass_inv[type] = 1.0/m;
  return true; //WARNING  
}

bool Atom_data::add_charges (unsigned int type, Real_t c) {
  if (type + 1 > owned.charge.size())
    owned.charge.resize (type + 1);
  owned.charge[type] = c;
  return true; //WARNING
}

void Atom_data::remove_atom(const int i) {
    owned.position.erase (owned.position.begin()+i);  
    owned.velocity.erase (owned.velocity.begin()+i);  
    owned.acceleration.erase (owned.acceleration.begin()+i);  
    owned.type.erase (owned.type.begin()+i);  
    owned.id.erase (owned.id.begin()+i);
    owned.msd_domain_cross.erase(owned.msd_domain_cross.begin()+i);
  --num_local_atoms;   
}

void Atom_data::remove_atom(std::vector<int> v_delete_list) {
  // sort them from greatest to lowest to maintain lower index
  std::sort(v_delete_list.begin(), v_delete_list.end(), std::greater<int>()); 
  for (auto i : v_delete_list) {
    owned.position.erase (owned.position.begin()+i);  
    owned.velocity.erase (owned.velocity.begin()+i);  
    owned.acceleration.erase (owned.acceleration.begin()+i);  
    owned.type.erase (owned.type.begin()+i);  
    owned.id.erase (owned.id.begin()+i);
    owned.msd_domain_cross.erase(owned.msd_domain_cross.begin()+i);
    --num_local_atoms;   
  }
}

void Atom_data::add_random_velocity() {

  // setting random velocity to all the particles
  std::mt19937 mt(1);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (auto && v : owned.velocity) {
    v.x = dist(mt);
    v.y = dist(mt);
    v.z = dist(mt);
  }


  // removing any center of mass velocity
  auto v_cm = owned_velocity_cm();
  auto psize = owned.velocity.size();
  for (unsigned int i = 0; i < psize; ++i) {
    if (owned.mass[owned.type[i]] != 0.0)
      owned.velocity[i] -= v_cm;
  }
}



void Atom_data::initialize_reading_xyz_frames(std::string input_file_name) {
  std::cout << "ai 1 : input_file_name : " << input_file_name << std::endl;
  ifs_xyz_postprocess.open(input_file_name.c_str());
  std::cout << "ai 2 : input_file_name : " << input_file_name << std::endl;
}
  
void Atom_data::finalize_reading_xyz_frames() {
  ifs_xyz_postprocess.close();
}   

int Atom_data::read_next_xyz_frame (bool set_frame, bool read_velocity) {
  

  auto & ifs = ifs_xyz_postprocess;



  int num_of_atoms = 0;
  
  ifs >> num_of_atoms;
  
  if (ifs.eof()) return -1; // don't repeat the last line

  {
    // Ignore the second line    
    std::string dummyLine;
    getline(ifs, dummyLine);
    getline(ifs, dummyLine);
  } 

  std::cout << "====================\n====================\n";
  
  for (int j = 0; j < num_of_atoms; ++j) {
  
    int type;
    
    ifs >> type;
    
    double x = 0, y = 0, z = 0, vx = 0, vy = 0, vz = 0;        
    
    ifs >> x >> y >> z;
    
    if (read_velocity)
        ifs >> vx >> vy >> vz;
    
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      
    std::cout <<type << " " << x << " " << y << " " << z <<  "\n";              
    
    if (set_frame)
    {
      owned.position[j].x = x;
      owned.position[j].y = y;
      owned.position[j].z = z;
      
      if (read_velocity) {
        owned.velocity[j].x = vx;
        owned.velocity[j].y = vy;
        owned.velocity[j].z = vz;
      }
    }
    
  }

    


  return 0;
}

} //objects

} // namespace caviar


