
//========================================================================
//
// Copyright (C) 2019 by deal.II authors and
// Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Some part of this file (Solving Laplace equation) is based on the step-40
// tutorial program of the deal.II library, with extensive modifications
// by Morad Biagooi and Ehsan Nedaaee Oskoee.
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_PLTDEALIIMPI_H
#define CAVIAR_OBJECTS_FORCEFIELD_PLTDEALIIMPI_H

#include "caviar/objects/force_field.h"

#ifdef CAVIAR_WITH_DEALII_MPI

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
//#include <deal.II/numerics/matrix_all_interpreter_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

// -- new inclusion.


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>



namespace LA {
#if defined(DEAL_II_WITH_PETSC) && !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>



CAVIAR_NAMESPACE_OPEN

class Atom_data;
namespace force_field {
using namespace dealii;

/**
 * This class calculates the effect of conductive boundaries on the charged
 * particles in MPI mode.
 */
class Plt_dealii_mpi : public Force_field
{
public:
  Plt_dealii_mpi (class CAVIAR *);
  ~Plt_dealii_mpi ();

  void calculate_acceleration ();  
  bool read (caviar::interpreter::Parser *);
  void verify_settings ();

  double total_potential_of_charges (const dealii::Point<3> &p);
  
public:
  void test_code();

  void run ();

  void read_domain();

  void make_boundary_face_normals ();

  // simple: a dealii standard Laplace solving process is done in every step.
  // faster: Laplace is solved using dealii::Filtered matrix, in order to get
  // rid of setup and assemble part of FE process in every step.
  // global: the mesh is not refined or refined only globally
  // adaptive: the mesh is 
  // hanging_nodes - partial_refine
  // filtered_matrix_hanging_nodes - partial_refine faster
  // simple_filtered_matrix _ simple and faster

  // simple_global : 
  // simple_adaptive :
  // faster_global :
  // faster_adaptive :

  void sa_setup_system();
  void sa_assemble_system();
  void sa_solve ();


  void calculate_induced_charge (const int);

  void output_vtk_solution (const int) const;

  void calculate_all_particles_mesh_force_acc();



  int tot_no_matched, tot_no_corrected;                                              
    
  parallel::distributed::Triangulation<3> triangulation;
  parallel::distributed::Triangulation<3> tria_reserve;


  IndexSet                                  locally_owned_dofs;
  IndexSet                                  locally_relevant_dofs;

  LA::MPI::SparseMatrix                     system_matrix;
  LA::MPI::Vector                           locally_relevant_solution;
  LA::MPI::Vector                           system_rhs;



  FE_Q<3> fe;
  DoFHandler<3> dof_handler;

  ConstraintMatrix     constraints;  

  int num_quadrature_points;
  double solver_control_tolerance;  
  bool use_preconditioner;

  class caviar::Atom_data * atom_data;

  bool initialized;
  
  std::vector<std::pair<int,double>> boundary_id_value;
  std::vector<unsigned> boundary_id_list;
  std::vector<int> refine_sequence_type, refine_sequence_value;  
  double k_electrostatic;
  std::vector<std::string> unv_mesh_filename;  
  int time_step_count, time_step_solve, time_step_output_vtk;
  int time_step_induced_charge;
  bool output_vtk;
  double derivation_length;

  std::vector<Tensor<1,3,double>> face_normal;
  std::vector<Point<3>> face_center;  
  std::vector<double> face_area;
  std::vector<unsigned> face_id;  
  std::vector<unsigned> face_id_ignore;
  std::ofstream ofs_induced_charge;
  bool induced_charge_id_init;
  bool output_induced_charge;  
  unsigned int boundary_id_max;
  bool make_time_profile;
  int zeros_of_mesh_output;
  int my_mpi_rank, mpi_world_size;

  std::vector<Vector<double>> face_center_pos, face_center_field;
  std::vector<double> face_center_potential;

  std::vector<caviar::Force_field*> force_field_custom;

  bool ignore_point_out_of_mesh;

};



//==================================================
//==================================================
//==================================================

namespace plt_dealii_mpi {

class BoundaryValues : public Function<3>
{
public:
  BoundaryValues () : Function<3>() {}
  
  BoundaryValues (double tp) : Function<3>(), total_potential{tp} {}  
  
  BoundaryValues (double tp,
    class caviar::force_field::Plt_dealii_mpi* df )
    : 
    Function<3>(),
    total_potential{tp},
    deal_force{df}  {

    }
       
  virtual double value (const Point<3> &p,
                        const unsigned int component = 0) const;                    
public:
  double total_potential;
  double potential_of_free_charges  (const dealii::Point<3> &p) const;

  class caviar::force_field::Plt_dealii_mpi *deal_force;
};
} //plt_dealii_mpi
} //force_field

CAVIAR_NAMESPACE_CLOSE

#else

CAVIAR_NAMESPACE_OPEN

namespace force_field {

class Plt_dealii_mpi : public Force_field
{
public:
  Plt_dealii_mpi (class CAVIAR *fptr);
  ~Plt_dealii_mpi ();
  void calculate_acceleration ();
  bool read (caviar::interpreter::Parser *);

};
} //finite_element

CAVIAR_NAMESPACE_CLOSE

#endif
#endif
