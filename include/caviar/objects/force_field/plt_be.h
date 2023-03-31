
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_PLTBE_H
#define CAVIAR_OBJECTS_FORCEFIELD_PLTBE_H

#include "caviar/objects/force_field.h"


namespace caviar {

class Atom_data;
namespace shape {
class Polyhedron;
}
namespace force_field {


/**
 * This class calculates the effect of conductive boundaries on the charged
 * particles. 
 */

class Plt_be : public Force_field
{
public:
  Plt_be (class CAVIAR *);
  ~Plt_be ();

  void calculate_acceleration ();  
  void verify_settings ();
  bool read (caviar::interpreter::Parser *);

  caviar::shape::Polyhedron *polyhedron;
  
  double k_electrostatic;

  std::vector<caviar::Force_field*> force_field_custom;

  //double total_potential_of_charges (const dealii::Point<3> &p);
  
  void start_spherical_test();
  void write_spherical_test();

  caviar::Force_field *test_force_spherical;
  std::ofstream ofs_test_force_spherical;
  bool init_test_force_spherical;
  std::string spherical_test_file_name;
  std::ofstream ofs_induced_charge;

  void jacobian_calculation();
  void make_inverse_matrix(); //  makes 'm_inverse' of 'D1'. Once per geometry.

  void set_potential_on_boundary();
  void make_vec_zz(); // makes 'vec_zz' from 'D2', 'm_inverse', and 'phi_boundary'.

  double potential_value(const Vector<double> v);

  std::vector<double> vec_zz;
  std::vector<double> phi_boundary;
  std::vector<double> jac;

  std::vector<std::vector <double> > D1, D2;

  std::vector<Vector<double>> pc1, pc2, pc3; // ?
  //std::vector<double> D_1, D_2;
  std::vector<Vector<double>> face_center ; // faces centeral point ?
  std::vector<double> tg, vg;
  unsigned face_size;

  //std::vector<std::pair<int,double>> boundary_id_value;
  std::vector<double> boundary_id_value;

  std::vector<unsigned> boundary_id_list;
  
  std::vector<std::vector<double>> m_inverse;
  bool initialized;

/*

  void run ();
  void run_time_profile ();
  void read_domain();
  void make_spherical_grid ();
  void make_boundary_face_normals ();
  void output_boundary_id_areas ();

  double calculate_induced_charge (const int, const int requested_id = -1);

  void output_vtk_solution (const int) const;
  void output_vtk_solution (const std::string filename) const;

  void calculate_all_particles_mesh_force_acc();

  void set_spherical_manifold();


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
  

  std::vector<int> refine_sequence_type, refine_sequence_value;  

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

  bool induced_charge_id_init;
  bool output_induced_charge;  
  unsigned int boundary_id_max;
  bool make_time_profile;
  int zeros_of_mesh_output;

  std::vector<Vector<double>> face_center_pos, face_center_field;
  std::vector<double> face_center_potential;



  bool ignore_point_out_of_mesh;
  bool time_profile;


  */
};


} //force_field

} // namespace caviar

#endif
