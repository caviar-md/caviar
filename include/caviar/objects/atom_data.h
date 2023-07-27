
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
#include "caviar/objects/atom_data/utility/bond.h"
#include "caviar/objects/atom_data/utility/angle.h"
#include "caviar/objects/atom_data/utility/proper_dihedral.h"

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
namespace neighborlist
{
  class Cell_list;
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

  /**
   * Initial setting of number of atom types.
   */
  virtual void set_num_atom_types(AtomType_t n) { num_atom_types = n; }

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
   * a simple adding a random velocity to all of the atoms. It may be useful
   * in some initial system settings.
   */
  virtual void add_random_velocity();

  /**
   * find and exchange owned atoms between domain or do periodic exchange
   */
  virtual bool exchange_owned();

  /**
   * find and exchange ghost atoms between domains or do periodic ghost
   */
  virtual void exchange_ghost();

  /**
   * does as it says
   */
  virtual bool position_inside_local_domain(const Vector<double> &pos);

  /**
   * finds the global_id of an Atom that's going to be added.
   */
  virtual GlobalID_t get_global_id();

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
   * Here we define an unnamed interpreter and make two of it. It can have a name
   * but since it is used only once here, its name won't be of any use. Also
   * a good name would be 'atom_data' which is used before!.
   */
  struct
  {
    /**
     * 'id' is a global and unique number assigned to an atom.
     */
    std::vector<GlobalID_t> id;

    /**
     * A tag to be done on the atoms.
     */
    std::vector<AtomType_t> tag;

    /**
     * Atom type decides the charge, mass and any other property shared between
     * a defined type (for example, Elements).
     */
    std::vector<AtomType_t> type;

    /**
     * 'mass' of an atom defined by the type. The mass may be used in
     * center-of-mass calculations and other functions. Do not depercate it.
     */
    std::vector<Real_t> mass;

    /**
     * simply the inverse value of 'mass' of an atom defined by the type.
     * since mass inverse is used in acceleration calculations.
     *
     */
    std::vector<Real_t> mass_inv;

    /**
     * 'charge' of an atom defined by the type.
     */
    std::vector<Real_t> charge;

    /**
     * 'radius' of an atom defined by the type. The user and the developers are
     * free to use this variable (for now!).
     */
    std::vector<Real_t> radius;

    /**
     * Different by atom_id. Can be changed in simulation (it is needed to be sent-recv. by MPI)
     */
    std::vector<Real_t> charge_atom;

    /**
     * Different by atom_id. Can be changed in simulation (it is needed to be sent-recv. by MPI)
     */
    std::vector<Real_t> mass_atom;

    /**
     * Atom kinematic properties in the current time-step.
     */
    std::vector<Vector<Real_t>> position, velocity, acceleration;

    /**
     * This vectors are used in some integrator schemes and constraint methods.
     * They can be defined in their related objects, but they may be needed in
     * more than one objects at once (for example, constraint::M_shake and
     * integrator::Leap_frog). This makes it the reason to define it here.
     * This function may be needed to have MPI_send-recv. process in these case.
     * look up to it.
     */
    std::vector<Vector<Real_t>> position_old, velocity_old, acceleration_old;

    /**
     * this vector is meaningful when there's one domain. We can calculate MSD
     * using this. It collects number of periodic domain cross for each particle.
     */
    std::vector<Vector<int>> msd_domain_cross;

    /**
     * this vector contain a molecule index for all the atoms. if it's '-1' the
     * atom is not of any molecule. This matters in the MPI process. All of the
     * atoms of a molecule should be existed in one process.
     */
    std::vector<int> molecule_index; //

    /**
     * number of total molecules.
     */
    int num_molecules;

    /**
     * The first std::vector, is the molecule index. the inner data contain bonds.
     */
    std::vector<std::vector<atom_data::Bond>> atomic_bond_vector;

    /**
     * Number of atomic bonds each atom have. It is used to limit
     * the bond creations.
     */
    std::vector<int> atomic_bond_count;

    /**
     * The first index, meaning
     * the first std::vector, is the molecule index. the inner data contain angles.
     */
    std::vector<std::vector<atom_data::Angle>> atomic_angle_vector;

    std::vector<std::vector<atom_data::Proper_dihedral>> atomic_properdihedral_vector;

  }
  /**
   * There are two types of this unnamed interpreter: owned and ghost. 'owned' atoms
   * are the one that matter in integrators in the domain. Ghost particles are
   * the particle near the domain boundaries which are not from the domain. Their
   * usage is for short-range force-field calculations.
   */
  owned,
      ghost;

  /**
   * it turns the process of recording owned data in the releated std::vector
   */
  bool record_owned_position_old, record_owned_velocity_old, record_owned_acceleration_old;

  /**
   * since not all of the situations need ghost particles velocity, we only send
   * if it is required.
   */
  bool make_ghost_velocity;

  /**
   * what these variables do are obvious. 'est' is for estimation.
   */
  LocalID_t num_local_atoms, num_local_atoms_est;
  GlobalID_t num_total_atoms;
  AtomType_t num_atom_types;

  std::vector<int> ghost_rank; // the rank of the domain in which the owned counterpart exists

  /**
   * i wonder if this should be here or it should be in the neighborlist class
   */
  std::vector<Vector<Real_t>> last_reneighborlist_pos;

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
   * Add mean square distance (MSD) calculations if needed. The default value is false for performance.
   */
  bool msd_process = false;

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

  /**
   * usage in 'empty_of_atoms()' functions.
   */
  class neighborlist::Cell_list *cell_list;

  FC_BASE_OBJECT_COMMON_TOOLS
};

CAVIAR_NAMESPACE_CLOSE

#endif
