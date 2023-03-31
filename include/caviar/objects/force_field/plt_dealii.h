
//========================================================================
//
// Copyright (C) 2019 by deal.II authors and
// Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Some part of this file (Solving Laplace equation) is based on the step-6
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_PLTDEALII_H
#define CAVIAR_OBJECTS_FORCEFIELD_PLTDEALII_H

#include "caviar/objects/force_field.h"

#ifdef CAVIAR_WITH_DEALII

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

namespace caviar {

class Atom_data;
namespace unique
{
class Time_function;
class Time_function_3d;
}
namespace force_field {
using namespace dealii;

/**
 * This class calculates the effect of conductive boundaries on the charged
 * particles. 
 */
class Plt_dealii : public Force_field
{
public:
  Plt_dealii (class CAVIAR *);
  ~Plt_dealii ();

  void calculate_acceleration ();  
  void verify_settings ();
  bool read (caviar::interpreter::Parser *);

  double total_potential_of_charges (const dealii::Point<3> &p);
  
public:



  void run ();
  void run_time_profile ();
  void read_domain();
  void make_spherical_grid ();
  void make_boundary_face_normals ();
  void output_boundary_id_areas ();
  void output_field_vectors(caviar::interpreter::Parser *);
  void output_potential_values(caviar::interpreter::Parser *);

  double potential (const Vector<double> &v); // Gives the total potential (sum of smooth and singular).
  double potential (const int);

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

  // simple_global; solve_type = 0;
  void sg_setup_system();
  void sg_assemble_system();
  void sg_solve ();

  // simple_adaptive; solve_type = 1;
  void sa_setup_system();
  void sa_assemble_system();
  void sa_solve ();

  // faster_global; solve_type = 2;
  void fg_setup_system();
  void fg_assemble_system();
  void fg_solve ();

  // faster_adaptive; solve_type = 3;
  void fa_setup_system();
  void fa_assemble_system();
  void fa_solve ();
  void fa_solve_time_profile ();


  double calculate_induced_charge (const int, const int requested_id = -1);

  void output_vtk_solution (const int) const;
  void output_vtk_solution (const std::string filename) const;

  void calculate_all_particles_mesh_force_acc();

  void set_spherical_manifold();

  void start_spherical_test();
  void write_spherical_test();

  void generate_ml_training_data(caviar::interpreter::Parser *parser);
 
  Triangulation<3>   tria_reserve;  

                                               
    
  Triangulation<3> triangulation;
  FE_Q<3> fe;
  DoFHandler<3> dof_handler;

  SparsityPattern sparsity_pattern;
  SparseMatrix<double> system_matrix;
  
  ConstraintMatrix     constraints;

  dealii::Vector<double> solution;
  dealii::Vector<double> system_rhs;


  int num_quadrature_points;
  int solver_control_maximum_iteration;  
  double solver_control_tolerance;  

  // solve_type related:
  int solve_type; // 0 1 2 3
  bool partial_refine_called; // use for hanging nodes warning (if boundary_refine or adaptive_refine called.
  bool re_do_setup_assemble; // used for 'faster' option, when there's the need to do setup and assemble part.
  bool use_preconditioner;
  double preconditioner_relaxation; // as is said in dealii docs,  0 < preconditioner_relaxation < 2
 
  class caviar::Atom_data * atom_data;

  bool initialized;
  
  unique::Time_function_3d *position_offset = nullptr;

  std::vector<std::pair<int,double>> boundary_id_value;
  std::vector<std::pair<int,caviar::unique::Time_function*>> boundary_id_time_function;
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
#if defined(CAVIAR_WITH_MPI)
  int my_mpi_rank, mpi_world_size;
#endif  

  std::vector<Vector<double>> face_center_pos, face_center_field;
  std::vector<double> face_center_potential;

  std::vector<caviar::Force_field*> force_field_custom;

  bool ignore_point_out_of_mesh;
  bool time_profile;

  caviar::Force_field *test_force_spherical;
  std::ofstream ofs_test_force_spherical;
  bool init_test_force_spherical;
  std::string spherical_test_file_name;

};



//==================================================
//==================================================
//==================================================

namespace plt_dealii {

class BoundaryValues : public Function<3>
{
public:
  BoundaryValues () : Function<3>() {}
  
  BoundaryValues (double tp) : Function<3>(), total_potential{tp} {}  
  
  BoundaryValues (double tp,
    class caviar::force_field::Plt_dealii* df )
    : 
    Function<3>(),
    total_potential{tp},
    deal_force{df}
    {}
       
  virtual double value (const Point<3> &p,
                        const unsigned int component = 0) const;                    
public:
  double total_potential;
  double potential_of_free_charges  (const dealii::Point<3> &p) const;

  class caviar::force_field::Plt_dealii *deal_force;

};
} //plt_dealii
} //force_field

} // namespace caviar

#else

namespace caviar {

namespace force_field {

class Plt_dealii : public Force_field
{
public:
  Plt_dealii (class CAVIAR *fptr);
  ~Plt_dealii ();
  void calculate_acceleration ();
  bool read (caviar::interpreter::Parser *);

};
} //finite_element

} // namespace caviar

#endif
#endif
