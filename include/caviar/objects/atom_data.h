
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

#ifndef CAVIAR_OBJECTS_ATOMDATA_H
#define CAVIAR_OBJECTS_ATOMDATA_H

#include "caviar/utility/objects_common_headers.h"
#include "caviar/objects/atom_data/utility/atom_struct.h"
#include "caviar/objects/atom_data/utility/atom_type_params.h"
#include "caviar/objects/atom_data/utility/molecule_struct.h"
#include "caviar/objects/atom_data/utility/molecule_type_params.h"
#include "caviar/objects/atom_data/utility/mpi_packet_info.h"
#include "caviar/objects/atom_data/utility/mpi_optimization.h"


CAVIAR_NAMESPACE_OPEN

class Domain;
namespace unique
{
  class Atom;
  class Atom_group;
  class Atom_list;
  class Molecule;
  class Molecule_group;
  class Molecule_list;
  class Time_function_3d;
}


/**
 * This class is the base class for all the atom_datas.
 * Atom_data contains all of the molecular and atomic data for a MD simulation.
 * It also handles data exchange between MPI domains.
 */
class Atom_data : public Pointers
{
public:
  Atom_data(class CAVIAR *);
  virtual ~Atom_data();
  virtual bool read(class caviar::interpreter::Parser *);

  /**
   * It represents the position of the origin of non_inertia Cartesian reference frame by a time function.
   * It will be used in ???.
   */
  // unique::Time_function_3d *position_offset = nullptr;


   /* 
   * More sharing results in less MPI overhead and probably faster MPI simulation. However, simulation will takes more RAM.
   * For large number of processors and particles, it must be tested.
   */
  MpiOptimization mpiOptimization = MpiOptimization::None;


  /**
   * It represents the velocity of origin of the non_inertia Cartesian reference frame by a time function.
   * It will be used in calculation of Temperature and Kinetic energy.
   */
  unique::Time_function_3d *velocity_offset = nullptr;

  /**
   * check if a position is empty of any atom. Usage in fixing number of atoms
   * or molarity when one wants to create atoms.
   */
  virtual bool empty_of_atoms(const Vector<Real_t>, double radius);

  /**
   * MPI case: Finds the correct MPI domain of the atoms
   */
  void set_atoms_mpi_rank();

  /**
   * MPI case: Finds the correct MPI domain of the atoms
   */
  int get_mpi_rank();


  /**
   * Import an xyz file contaning atoms positions and maybe velocities
   */
  bool add_xyz_data_file(caviar::interpreter::Parser *parser);

  /**
   * checks by atom_type
   */
  virtual bool empty_of_atoms(const Vector<Real_t>, int type);
  virtual bool empty_of_atoms(unique::Atom &a);
  virtual bool empty_of_atoms(unique::Molecule &m);

  /**
   *  position of the center of mass
   */
  virtual Vector<Real_t> owned_position_cm();

  /**
   *  velocity of the center of mass
   */
  virtual Vector<Real_t> owned_velocity_cm();

  /**
   *  angular momentum of the center of mass
   */
  virtual Vector<double> owned_angular_momentum_cm()
  {
    return owned_angular_momentum_cm(owned_position_cm());
  }

  /**
   *  angular momentum of the center of mass
   */
  virtual Vector<double> owned_angular_momentum_cm(const Vector<double> &p_cm);

  /**
   *  inertia_tensor of the center of mass. The tensor type may be modified
   *  in the future.
   */
  virtual std::vector<std::vector<double>> owned_inertia_tensor_cm()
  {
    return owned_inertia_tensor_cm(owned_position_cm());
  }

  /**
   *  inertia_tensor of the center of mass. The tensor type may be modified
   *  in the future.
   */
  virtual std::vector<std::vector<double>> owned_inertia_tensor_cm(const Vector<double> &p_cm);

  /**
   * Initial setting of number of atoms.
   */
  virtual void set_num_total_atoms(GlobalID_t);

  // /**
  //  * Initial setting of number of atom types.
  //  */
  // virtual void set_num_atom_types(AtomType_t n) { num_atom_types = n; }

  /**
   * total number of system degree of freedom. For simple atomic simulations,
   * it returns '3*num_total_atoms'. For molecular simulations,
   * '3*num_total_atoms - atomic_bonds - atomic_angles' is returned.
   * It should be developed for special cases.
   */
  virtual int degree_of_freedoms();

  /**
   * reserve the owned std::vector for a faster push_back assignment
   */
  virtual void reserve_owned_vectors();

  /**
   * gets the data and add it to the owned if it should be owned.
   */
  virtual bool add_atom(GlobalID_t,
                        AtomType_t,
                        const Vector<Real_t> &,
                        const Vector<Real_t> &vel = Vector<Real_t>{0.0, 0.0, 0.0});

  /**
   * add unique::Atom to the owned data
   */
  virtual bool add_atom(caviar::unique::Atom &a);
  virtual bool add_atom(caviar::unique::Atom_group &a);
  virtual bool add_atom(caviar::unique::Atom_list &a);

  /**
   * add unique::Molecule to the owned data
   */
  virtual bool add_molecule(caviar::unique::Molecule &m);
  virtual bool add_molecule(caviar::unique::Molecule_group &m);
  virtual bool add_molecule(caviar::unique::Molecule_list &m);

  /**
   * merging two molecules by their molecule index
   */
  void merge_molecules(int molecule_index_1, int molecule_index_2);

  /**
   * adds a new bond between existing atoms in the  atom_data and merge molecules if possible
   */
  void add_atomic_bond(const atom_data::Bond &bond);

  /**
   * adds a new bond between existing atoms in the  atom_data and merge molecules if possible
   */
  void add_atomic_angle(const atom_data::Angle &angle);

  /**
   * adds a new bond between existing atoms in the  atom_data and merge molecules if possible
   */
  void add_atomic_properdihedral(const atom_data::Proper_dihedral &proper_dihedral);

  /**
   * remove atomic bond if it exist. Also remove atomic angles and proper dihedrals if the bond is used in them.
   */
  void remove_atomic_bond(const atom_data::Bond &bond);

  /**
   * remove atomic angle if it exist
   */
  void remove_atomic_angle(const atom_data::Angle &angle);

  /**
   * remove atomic angle if it exist
   */
  void remove_atomic_properdihedral(const atom_data::Proper_dihedral &proper_dihedral);

  /**
   * remove atomic bond if it exist. Also remove atomic angles and proper dihedrals if the bond is used in them.
   */
  bool check_atomic_bond_exist(const atom_data::Bond &bond);

  /**
   * remove atomic angle if it exist
   */
  bool check_atomic_angle_exist(const atom_data::Angle &angle);

  /**
   * remove atomic angle if it exist
   */
  bool check_atomic_properdihedral_exist(const atom_data::Proper_dihedral &proper_dihedral);

  /**
   * sets the mass of an atom type
   */
  virtual bool add_masses(unsigned int, Real_t);

  /**
   * sets the charge of an atom type
   */
  virtual bool add_charges(unsigned int, Real_t);

  /**
   * does as it says.
   */
  virtual void remove_atom(const int index);
  virtual void remove_atom(std::vector<int> index_list);

  /**
   * calculates the instantaneous temperature of all of the owned atoms.
   * by means of equipartition theorem, 'k = 1/2 * k_b * N_df * T'.
   * It can be defined differently, or for a subset of degree of freedoms.
   */
  virtual double temperature();

  /**
   * calculates total kinetic energy of the atoms
   */
  virtual double kinetic_energy();

  /**
   * calculates total kinetic energy of a type of atoms
   */
  virtual double kinetic_energy(const int);

  /**
   * a simple adding a random velocity to all of the atoms. The argument is the seed for random initialization.
   */
  virtual void add_random_velocity(unsigned int seed);

  /**
   * find and exchange owned atoms between domain or do periodic boundary condition movement.
   */
  virtual bool exchange_owned(long step = -1);



  /**
   * find and exchange ghost atoms between domains or do periodic boundary condition for ghosts.
   */

  virtual void exchange_ghost(long step = -1);

  /**
   * does as it says
   */
  virtual bool position_inside_local_domain(const Vector<double> &pos);

  /**
   * in MPI, it is the sum of local atoms.
   */
  virtual long get_num_of_atoms_local();

  /**
   * in MPI, it is the sum of local atoms.
   */
  virtual long get_num_of_atoms_global();

  /**
   * used when there are different MPI domains and all of them should know all
   * the atom_data
   */
  virtual void synch_owned_data(int type);

  /**
   * recording owned.position, velocity... in the owned.position_old ... if
   * necessary.
   */
  virtual void record_owned_old_data();

  /**
   *  it resets all the acceleration to a zero vector.
   */
  virtual void reset_owned_acceleration();

  /**
   * This function is called before reading an xyz file,
   * frame by frame. It is developed for postprocessing an
   * xyz file.
   */
  virtual void initialize_reading_xyz_frames(std::string input_file_name);

  /**
   * This function is called after reading an xyz file,
   *
   *
   */
  virtual void finalize_reading_xyz_frames();

  /**
   * used by 'read_next_xyz_frame'. It is initialized in 'initialize_reading_xyz_frames'
   *
   */
  std::ifstream ifs_xyz_postprocess;

  /**
   * This function reads the next frame of xyz file.
   * it only sets the frame into atom_data if the 'set_frame'
   * argument is true.
   */
  virtual int read_next_xyz_frame(bool set_frame, bool read_velocity);

 
  /**
   *  reserving memory for atom_struct_owned to have faster atom import or atom addition
   */
  virtual bool atom_struct_owned_reserve(long new_size);

   /**
   *  resizing atom_struct_owned correctly. It will resize only the important params according to atom_data configuration.
   */
  virtual bool atom_struct_owned_resize(long new_size);

   /**
   *  resizing atom_struct_owned correctly. It will resize only the important params according to atom_data configuration.
   */
  virtual bool atom_struct_ghost_resize(long new_size);
  /**
   * 'owned' atoms are the one that matter in integrators in the domain. 
   */
  atom_data::Atom_struct atom_struct_owned;

  /**
   *  Ghost particles are the particle near the domain boundaries which are not from the domain. Their
   * usage is for short-range force-field calculations.
   */
  atom_data::Atom_struct atom_struct_ghost;

  /**
   * Atom type physical parameters. 
   */
  atom_data::Atom_type_params atom_type_params;

  /**
   * 'owned' atoms are the one that matter in integrators in the domain. 
   */
  std::vector<atom_data::Molecule_struct> molecule_struct_owned;

  
  /**
   * number of total molecules.
   */
  int num_molecules;

  /**
   * add_atom_id_to_molecule if it is not found.
   */
  void add_atom_id_to_molecule(int atom_id, int molecule_index);

  /**
   * remove_atom_id_to_molecule if it exist. if it doesn't exist, call error.
   */
  void remove_atom_id_from_molecule(int atom_id, int molecule_index);

  /**
   * index: the position of atom in std::vector.
   * id: unique identifier of atoms between all MPI processes.
   * In order of less MPI communications for molecules, this vector is implemented.
   * the id is global but index is local in MPI domain.
   * in serial mode, 'atom_id_to_index[id] = id' is possible.
   */
  std::vector<int> atom_id_to_index;

  /**
   * since not all of the situations need ghost particles velocity, we only send
   * if it is required.
   */
  bool make_ghost_velocity;

  // /**
  //  * what these variables do are obvious. 'est' is for estimation.
  //  */
  // LocalID_t num_local_atoms, num_local_atoms_est;
  // GlobalID_t num_total_atoms;
  // AtomType_t num_atom_types;


  /**
   * if true, more than just atom position have to be synched in single domain mpi case
   */
  bool synch_owned_data_bcast_details;

  /**
   * is it useful?
   */
  double neighborlist_cutoff;

  /**
   * obvious
   */
  double ghost_cutoff;

  /**
   * is it useful?
   */
  double cutoff_extra;

  /**
   * Boltzman constant. Used in temperature calculation.
   */
  double k_b;



  /**
   *  number of external degrees of freedom according to page 114 of
   *  Philippe H. Hunenberger, Adv. Polym. Sci. (2005) 173:105–149  ,
   *  'N_r = 0 in the presence of stochastic and frictional forces.
   *  N_r = 3 under periodic boundary conditions,
   *  N_r = 6 under vacuum boundary conditions'
   */
  int n_r_df;

  /**
   *  If n_r_df is not set, it will be calculated and returned.
   */
  int get_n_r_df();



  /**
   * stochastic and frictional forces presence. It affects n_r_df.
   * If it is not activated, the system is under vacuum boundary condition or
   * periodic boundary condition.
   * This will be checked if 'n_r_df' has not been set.
   */
  bool stochastic_force_present;

  /**
   * usage
   */
  class Domain *domain;

  bool get_msd_process()
  {
    return msd_process;
  }

  void set_msd_process(bool stat)
  {
    msd_process = stat;
    atom_struct_owned.msd_domain_cross.resize(get_num_of_atoms_local(),caviar::Vector<int>{0,0,0});  
  }

  void set_record_owned_position_old(bool stat)
  {
    record_owned_position_old = stat;
    atom_struct_owned.position_old.resize(get_num_of_atoms_local(),caviar::Vector<double>{0,0,0});  
  }

    void set_record_owned_velocity_old(bool stat)
  {
    record_owned_velocity_old = stat;
    atom_struct_owned.velocity_old.resize(get_num_of_atoms_local(),caviar::Vector<double>{0,0,0});  
  }

    void set_record_owned_acceleration_old(bool stat)
  {
    record_owned_acceleration_old = stat;
    atom_struct_owned.acceleration_old.resize(get_num_of_atoms_local(),caviar::Vector<double>{0,0,0});  
  }
  //===========================
  // Private Section
  //===========================
  private:
  /**
   * The tools meant to use in MPI mode. In order to reduce memory allocation in each timestep, we define them out of scope
   */
  struct
  {
    std::vector<int> send_index[3][3][3]; // the index of std::vector<> of the owned
    std::vector<int> send_index_all;      // the index of std::vector<> of the owned, in a 1D vector

    std::vector<double> send_data[3][3][3]; 
    std::vector<double> recv_data[3][3][3];

    int send_num[3][3][3]; // num of owned to be send to the domain all[i][j][k]
    int recv_num[3][3][3]; // num of owned to be recieved from domain all[i][j][k]/ used in he

    int send_mpi_tag[3][3][3]; // since there might be two messages from the same domain to another but from different angles,
    int recv_mpi_tag[3][3][3]; // , this tag helps to distinguish messages form each other.

    bool initialize = true;

  } mpi_tools;

  /**
   * Add mean square distance (MSD) calculations if needed. The default value is false for performance.
   */
  bool msd_process = false;

  /**
   * it turns the process of recording owned data in the releated std::vector
   */
  bool record_owned_position_old, record_owned_velocity_old, record_owned_acceleration_old;

  /**
   * It stores the info about location of the owned Atom_struct data inside the MPI packets.
   */
  atom_data::MPI_packet_info mpi_packet_info_owned;

  /**
   * It stores the info about location of the ghost Atom_struct data inside the MPI packets.
   */
  atom_data::MPI_packet_info mpi_packet_info_ghost;
  
  /**
   * MPI rank of the classs
  */
  int my_mpi_rank = -1;

    /**
   * find and exchange owned atoms between domain or do periodic boundary condition movement.
   * It is used in serial simulation case or in MPI when we consider all MD calculations 
   * are in mpi_rank==0
   */
  virtual bool exchange_owned_single_md_domain(long step = -1);


  /**
   * find and exchange owned atoms between domain or do periodic boundary condition movement.
   * Memory efficient case.
   */
  virtual bool exchange_owned_mpi(long step = -1);

  /**
   * find and exchange owned atoms between domain or do periodic boundary condition movement.
   * CPU efficient case.
   */
  virtual bool exchange_owned_mpi_shared_atoms(long step = -1);

  /**
   * find and exchange ghost atoms between domain or do periodic boundary condition movement.
   * It is used in serial simulation case or in MPI when we consider all MD calculations 
   * are in mpi_rank==0
   */
  virtual void exchange_ghost_single_md_domain(long step = -1);


  /**
   * find and exchange ghost atoms between domain or do periodic boundary condition movement.
   * Memory efficient case.
   */
  virtual void exchange_ghost_mpi(long step = -1);

  /**
   * find and exchange ghost atoms between domain or do periodic boundary condition movement.
   * CPU efficient case.
   */
  virtual void exchange_ghost_mpi_shared_atoms(long step = -1);
  //===========================
  // Public Section
  //===========================
  public:

  FC_BASE_OBJECT_COMMON_TOOLS
};

CAVIAR_NAMESPACE_CLOSE

#endif
