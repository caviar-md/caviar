
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
#include "caviar/objects/unique.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/neighborlist/cell_list.h"
#include "caviar/objects/neighborlist/verlet_list.h"
#include "caviar/objects/unique/atom.h"
#include "caviar/objects/unique/atom_group.h"
#include "caviar/objects/unique/atom_list.h"
#include "caviar/objects/unique/molecule.h"
#include "caviar/objects/unique/molecule_group.h"
#include "caviar/objects/unique/molecule_list.h"

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

//void Atom_data::allocate () {
  
//}

//void Atom_data::output_data (int ) {
  
//}


/*
#define FIND_UNIQUE_OBJECT_BY_NAME_NIC(OBJECT_TYPE,ITERATOR_NAME) \
  std::string name_to_find_;\
  MAKE_ERROR_MASSAGE_EXPECTED(err,OBJECT_TYPE,"name.")\
  GET_A_STRING(name_to_find_,"", err)\
  CHECK_NAME_EXISTANCE(name_to_find_, ITERATOR_NAME, "","")\
  if (ITERATOR_NAME->second.type != caviar::interpreter::object_handler::gdst( #OBJECT_TYPE ))\
    error->all(FC_FILE_LINE_FUNC_PARSE,": undefined object. ");

#define FIND_UNIQUE_OBJECT_BY_NAME(OBJECT_TYPE,ITERATOR_NAME) \
  std::map<std::string,caviar::interpreter::object_handler::Dictionary>::iterator ITERATOR_NAME;\
  FIND_UNIQUE_OBJECT_BY_NAME_NIC(OBJECT_TYPE,ITERATOR_NAME)
*/
bool Atom_data::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"ghost_cutoff")) {
      GET_OR_CHOOSE_A_REAL(ghost_cutoff,"","")
      if (ghost_cutoff < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "ghost_cutoff have to non-negative."); 
    } else if (string_cmp(t,"neighborlist_cutoff")) {
      GET_OR_CHOOSE_A_REAL(neighborlist_cutoff,"","")
      if (neighborlist_cutoff < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "neighborlist_cutoff have to non-negative."); 
    } else if (string_cmp(t,"cutoff_extra")) {
      GET_OR_CHOOSE_A_REAL(cutoff_extra,"","")
      if (cutoff_extra < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "neighborlist_cutoff have to non-negative."); 
    } else if (string_cmp(t,"k_b")) {
      GET_OR_CHOOSE_A_REAL(k_b,"","")
      if (k_b < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "k_b have to non-negative."); 
    } else if (string_cmp(t,"add_atom")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,atom)
      objects::unique::Atom *a = dynamic_cast<objects::unique::Atom *>(object_container->unique[it->second.index]);
      add_atom(*a);
    } else if (string_cmp(t,"add_atom_group")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,atom_group)
      objects::unique::Atom_group *a = dynamic_cast<objects::unique::Atom_group *>(object_container->unique[it->second.index]);
      add_atom(*a);
    } else if (string_cmp(t,"add_atom_list")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,atom_list)
      objects::unique::Atom_list *a = dynamic_cast<objects::unique::Atom_list *>(object_container->unique[it->second.index]);
      add_atom(*a);
    } else if (string_cmp(t,"add_molecule")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,molecule)
      objects::unique::Molecule *a = dynamic_cast<objects::unique::Molecule *>(object_container->unique[it->second.index]);
      add_molecule(*a);
    } else if (string_cmp(t,"add_molecule_group")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,molecule_group)
      objects::unique::Molecule_group *a = dynamic_cast<objects::unique::Molecule_group *>(object_container->unique[it->second.index]);
      add_molecule(*a);
    } else if (string_cmp(t,"add_molecule_list")) {
      FIND_OBJECT_BY_NAME(unique,it)
      FC_CHECK_OBJECT_CLASS_NAME(unique,it,molecule_list)
      objects::unique::Molecule_list *a = dynamic_cast<objects::unique::Molecule_list *>(object_container->unique[it->second.index]);
      add_molecule(*a);
    } else if (string_cmp(t,"add_type_radius")) {
      //auto ind = parser->get_positive_int();
      //auto r = parser->get_real();
      int ind = 0;
      double r = 0;
      GET_OR_CHOOSE_A_INT(ind,"","")
      GET_OR_CHOOSE_A_REAL(r,"","")
      if (static_cast<int> (owned.radius.size()) < ind+1) {
        owned.radius.resize(ind+1);
      }
      owned.radius[ind] = r;
      if (r<0) output->warning("you have entered a negative value for atom radius.");
    } else if (string_cmp(t,"add_type_charge")) {
      //auto ind = parser->get_positive_int();
      //auto c = parser->get_real();
      int ind = 0;
      double c = 0;
      GET_OR_CHOOSE_A_INT(ind,"","")
      GET_OR_CHOOSE_A_REAL(c,"","")
      if (static_cast<int> (owned.charge.size()) < ind+1) {
        owned.charge.resize(ind+1);
      }
      owned.charge[ind] = c;
    } else if (string_cmp(t,"add_type_mass")) {

      //auto ind = parser->get_positive_int();
      //auto m = parser->get_real();
      int ind = 0;
      double m = 0;
      GET_OR_CHOOSE_A_INT(ind,"","")
      GET_OR_CHOOSE_A_REAL(m,"","")
      if (static_cast<int> (owned.mass.size()) < ind+1) {
        owned.mass.resize(ind+1);
        owned.mass_inv.resize(ind+1);
      }
      owned.mass[ind] = m;
      if (m==0.0) owned.mass_inv[ind] = 0.0;
      else owned.mass_inv[ind] = 1.0/m;
      if (m<0) output->warning("you have entered a negative value for atom mass.");
    } else if (string_cmp(t,"set_domain") || string_cmp(t,"domain")) {
      FIND_OBJECT_BY_NAME(domain,it)
      domain = object_container->domain[it->second.index];
    } else if (string_cmp(t,"set_cell_list") || string_cmp(t,"cell_list")) {
      FIND_OBJECT_BY_NAME(neighborlist,it)
      FC_CHECK_OBJECT_CLASS_NAME(neighborlist,it,cell_list)
      cell_list = dynamic_cast<objects::neighborlist::Cell_list *>(object_container->neighborlist[it->second.index]);
    } else if (string_cmp(t,"add_xyz_data_file")) {
      add_xyz_data_file(parser);
      return true;
    } else if (string_cmp(t,"set_owned_position")) {
      auto ind = parser->get_int();
      auto x = parser->get_real();
      auto y = parser->get_real();
      auto z = parser->get_real();
      owned.position[ind] = Vector<Real_t> {x, y, z};
    } else if (string_cmp(t,"set_owned_velocity")) {
      auto ind = parser->get_int();
      auto x = parser->get_real();
      auto y = parser->get_real();
      auto z = parser->get_real();
      owned.velocity[ind] = Vector<Real_t> {x, y, z};
    } else if (string_cmp(t,"set_owned_acceleration")) {
      auto ind = parser->get_int();
      auto x = parser->get_real();
      auto y = parser->get_real();
      auto z = parser->get_real();
      owned.acceleration[ind] = Vector<Real_t> {x, y, z};
    } else if (string_cmp(t,"add_random_velocity")) {
      add_random_velocity();
      return true;
    } else if (string_cmp(t,"n_r_df")) {
      GET_OR_CHOOSE_A_INT(n_r_df,"","")
    } else FC_ERR_UNDEFINED_VAR(t)

  }
  return in_file;
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

  auto sum_of_dihedrals = 0;
  for (auto && i : owned.atomic_properdihedral_vector)
    sum_of_dihedrals += i.size();


  return df - sum_of_bonds - sum_of_angles - sum_of_dihedrals - get_n_r_df();
}


Vector<Real_t> Atom_data::owned_position_cm () {
  Vector<Real_t> p_cm {0.0,0.0,0.0};
  double mass_sum = 0.0;
  auto p_size = owned.position.size(); // MPI check
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:p_cm,mass_sum)
#endif     
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
  auto p_size = owned.velocity.size();// MPI check
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:v_cm,mass_sum)
#endif     
  for (unsigned int i = 0; i < p_size; ++i) {
    auto type_i = owned.type[i];
    auto mass_i = owned.mass[type_i];
    mass_sum += mass_i;
    v_cm += owned.velocity[i] * mass_i;
  }
  v_cm = v_cm / mass_sum;
  return v_cm;
}

Vector<double> Atom_data::owned_angular_momentum_cm (const Vector<double> &p_cm) {
  
  Vector<double> L_cm {0.0, 0.0,0.0};
  
  auto p_size = owned.position.size();// MPI check
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:L_cm)
#endif     
  for (unsigned int i = 0; i < p_size; ++i) {
    auto type_i = owned.type[i];
    auto mass_i = owned.mass[type_i];
    
    L_cm += mass_i*cross_product(owned.position[i] - p_cm, owned.velocity[i]);    
  }  
  
  return L_cm;
}

std::vector<std::vector<double> > Atom_data::owned_inertia_tensor_cm (const Vector<double> &p_cm) {
  
  std::vector<std::vector<double> > I_cm (3, std::vector<double> (3, 0.0));
  
  auto p_size = owned.position.size(); // MPI check
#ifdef CAVIAR_WITH_OPENMP  
  #pragma omp parallel for reduction (+:it_cm)
#endif     
  for (unsigned int i = 0; i < p_size; ++i) {
    auto type_i = owned.type[i];
    auto mass_i = owned.mass[type_i];
    auto p = owned.position[i] - p_cm; //relative position
    I_cm[0][0] += mass_i*(p.y*p.y + p.z*p.z);
    I_cm[1][1] += mass_i*(p.x*p.x + p.z*p.z);    
    I_cm[2][2] += mass_i*(p.x*p.x + p.y*p.y);
    
    I_cm[0][1] -= mass_i*(p.x*p.y);
    I_cm[1][2] -= mass_i*(p.y*p.z);    
    I_cm[2][0] -= mass_i*(p.z*p.x);
    
    I_cm[1][0] = I_cm[0][1];
    I_cm[2][1] = I_cm[1][2];
    I_cm[0][2] = I_cm[2][0];
  }      
  
  return I_cm;
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


