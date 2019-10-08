
//========================================================================
//
// Copyright (C) 2019 by deal.II authors and
// Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Some part of this file (Solving Laplace equation) is based on the step-6
// tutorial program of the deal.II library, with extensive modifications
// by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// The CAVIAR package is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the CAVIAR distribution.
//
//========================================================================

#include "caviar/objects/force_field/plt_dealii.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/interpreter/communicator.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/utility/macro_constants.h"

#ifdef CAVIAR_WITH_DEALII

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

//#include <deal.II/grid/grid_all_interpreter_tools.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

//#include <deal.II/dofs/dof_all_interpreter_tools.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

//#include <deal.II/numerics/vector_all_interpreter_tools.h>
#include <deal.II/numerics/vector_tools.h>

//#include <deal.II/numerics/matrix_all_interpreter_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/filtered_matrix.h>

#ifdef CAVIAR_WITH_DEALII_WITH_OPENCASCADE
 #include <deal.II/opencascade/boundary_lib.h>
 #include <deal.II/opencascade/utilities.h>
#endif

#include <fstream>
#include <iostream>
#include <cmath>
#include <time.h>

#include "caviar/objects/force_field/utility/plt_dealii_constants.h"
#include "caviar/objects/force_field/utility/plt_dealii_functions.h"
   
namespace caviar {
namespace objects {
namespace force_field {

//==================================================
//==================================================
//==================================================



Plt_dealii::Plt_dealii(CAVIAR * fptr) : 
    caviar::objects::Force_field{fptr},
    fe (FC_DEALII_FE_Q_POLY_DEGREE),
    dof_handler (triangulation), 
    num_quadrature_points{FC_DEALII_FE_Q_POLY_DEGREE+1}, 
    solver_control_maximum_iteration{1000},
    solver_control_tolerance{1e-6}, 
    solve_type{3},
    partial_refine_called{false},
    re_do_setup_assemble{false},
    use_preconditioner{true},
    preconditioner_relaxation{1.2},
    atom_data{nullptr}, 
    initialized{false}, k_electrostatic{1.0},
    time_step_count{0}, time_step_solve{1},
    time_step_output_vtk{1}, time_step_induced_charge{1},
    output_vtk{false}, 
    induced_charge_id_init{false},
    output_induced_charge{false}, boundary_id_max{0},
    make_time_profile{false},
    zeros_of_mesh_output{2}
#if defined(CAVIAR_WITH_MPI)
    ,my_mpi_rank{fptr->comm->me}
    ,mpi_world_size{fptr->comm->nprocs}
#endif  
{ 
  FC_OBJECT_INITIALIZE_INFO
  ignore_point_out_of_mesh = false;
  time_profile = false;
  output->warning( " "
      "Always add force_field::DealII_poisson_custom after 'ewald_k' and"
      " 'slab' to the integrators because they have to be initialized in"
      " every steps. The initialization functions are done in "
      "calculate_acceleration function.");
  test_force_spherical=nullptr;
  init_test_force_spherical=false;
}

//==================================================
//==================================================
//==================================================

Plt_dealii::~Plt_dealii() {}

//==================================================
//==================================================
//==================================================


bool Plt_dealii::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"make_grid")) {
      //make_grid();
        error->all (FC_FILE_LINE_FUNC_PARSE, "not implemented");
    } else if (string_cmp(t,"num_quadrature_points")) {    
      num_quadrature_points = parser->get_literal_int();     
    } else if (string_cmp(t,"refine_global")) {    
      int n = parser->get_literal_int();     
      refine_sequence_type.push_back (0);
      refine_sequence_value.push_back (n);
    } else if (string_cmp(t,"refine_boundary")) {
      int n = parser->get_literal_int();    
      refine_sequence_type.push_back (1);
      refine_sequence_value.push_back (n);
    } else if (string_cmp(t,"boundary_id_value")) {
      int id = parser->get_literal_int();
      if (id == 0) {
        error->all (FC_FILE_LINE_FUNC_PARSE, "#boundary_id=0 is reserved for free boundary");
      }
      double value = 0;
      GET_OR_CHOOSE_A_REAL(value,"","")
      boundary_id_value.push_back(std::make_pair(id, value));
      boundary_id_list.push_back(id);
    }  else if (string_cmp(t,"k_electrostatic")) {
      GET_OR_CHOOSE_A_REAL(k_electrostatic,"","")    
      if (k_electrostatic < 0)  error->all (FC_FILE_LINE_FUNC_PARSE, "k_electrostatic has to be non-negative.");            
    } else if (string_cmp(t,"preconditioner_relaxation")) {
      preconditioner_relaxation = parser->get_literal_real();
      if (preconditioner_relaxation<0 || preconditioner_relaxation >2) 
        error->all (FC_FILE_LINE_FUNC, "expected '0.0 < preconditioner_relaxation < 2.0'.");
    } else if (string_cmp(t,"solver_control_tolerance")) {
      solver_control_tolerance = parser->get_literal_real();
    } else if (string_cmp(t,"add_unv_mesh")) {
      auto token = parser->get_val_token();
      auto file_name = token.string_value;
      unv_mesh_filename.push_back(file_name);
    } else if (string_cmp(t,"read_unv_mesh")) {
      dealii_functions::read_domain(triangulation, tria_reserve, unv_mesh_filename);
    } else if (string_cmp(t,"time_step_solve")) {
      time_step_solve = parser->get_literal_int();    
    } else if (string_cmp(t,"output_vtk")) {
      output_vtk = true;
      time_step_output_vtk = parser->get_literal_int();
    } else if (string_cmp(t,"output_induced_charge")) {
      output_induced_charge = true;
      time_step_induced_charge = parser->get_literal_int();
#if defined(CAVIAR_WITH_MPI)
      if (my_mpi_rank == 0)     
        ofs_induced_charge.open ("o_induced_charge");
#else 
      ofs_induced_charge.open ("o_induced_charge");
#endif
    } else if (string_cmp(t,"output_induced_charge_name")) {
      output_induced_charge = true;
      time_step_induced_charge = parser->get_literal_int();
      auto token = parser->get_val_token();
      std::string fn = token.string_value;
#if defined(CAVIAR_WITH_MPI)
      if (my_mpi_rank == 0)     
        ofs_induced_charge.open (fn.c_str());
#else 
      ofs_induced_charge.open (fn.c_str());
#endif
    } else if (string_cmp(t,"induced_charge_ignore_id")) {
      unsigned id = parser->get_literal_int();
      face_id_ignore.push_back(id);
    } else if (string_cmp(t,"dealii_grid_generator")) {
      dealii_functions::dealii_grid_generator(fptr, parser, triangulation);
      return in_file;
    } else if (string_cmp(t,"add_force_field") || string_cmp(t,"force_field")) {
      FIND_OBJECT_BY_NAME(force_field,it)
      force_field_custom.push_back(object_container->force_field[it->second.index]);
    } else if (string_cmp(t,"set_solve_type")) {
      auto t1 = parser->get_val_token();
      auto st = t1.string_value;
      if (string_cmp(st,"simple_global")) {
        solve_type = 0;
      } else if (string_cmp(st,"simple_adaptive")) {
        solve_type = 1;
      } else if (string_cmp(st,"faster_global")) {
        solve_type = 2;
      } else if (string_cmp(st,"faster_adaptive")) {
        solve_type = 3;
      } else {
        error->all(FC_FILE_LINE_FUNC, static_cast<std::string>("undefined option : ")+st);
      } 
    } else if (string_cmp(t,"re_initialize")) {
      initialized = false;
    } else if (string_cmp(t,"calculate_acceleration")) {
      calculate_acceleration();
    } else if (string_cmp(t,"ignore_point_out_of_mesh")) {
      ignore_point_out_of_mesh = true;
    } else if (string_cmp(t,"time_profile")) {
      time_profile = true;
    } else if (string_cmp(t,"write_spherical_test")) {
      write_spherical_test();
    } else if (string_cmp(t,"set_spherical_test_force")) {
      FIND_OBJECT_BY_NAME(force_field,it)
      test_force_spherical = object_container->force_field[it->second.index];
    } else if (string_cmp(t,"start_spherical_test")) {
      auto token = parser->get_val_token();
      spherical_test_file_name = token.string_value;
      ofs_test_force_spherical.open(spherical_test_file_name.c_str());
      start_spherical_test();

    } else if (string_cmp(t,"finish_spherical_test")) {
      ofs_test_force_spherical.close();
      ofs_induced_charge.close();
    } else if (string_cmp(t,"refine_global_here")) {
      triangulation.refine_global(1);
    } else if (string_cmp(t,"output_vtk_here")) {
      auto token = parser->get_val_token();
      spherical_test_file_name = token.string_value;
      output_vtk_solution (spherical_test_file_name);
    }
    
    else FC_ERR_UNDEFINED_VAR(t)
  }
  
  return in_file;
}
//==================================================
//==================================================
//==================================================

void Plt_dealii::start_spherical_test() {
  FC_NULLPTR_CHECK(test_force_spherical)
  FC_NULLPTR_CHECK(atom_data)
  if (force_field_custom.size()<1) 
    error->all(FC_FILE_LINE_FUNC,"expected an electrostatic forcefield added");
  initialized = false;
  //if (!init_test_force_spherical) {
  //  init_test_force_spherical = true;
  //}
  make_boundary_face_normals (); 
  output_boundary_id_areas ();

    switch(solve_type) {
    case 0 :
      sg_setup_system ();
      sg_assemble_system ();
      sg_solve ();      
    break;

    case 1 :
      sa_setup_system ();
      sa_assemble_system ();
      sa_solve ();      
    break;

    case 2 :
      fg_setup_system ();
      fg_assemble_system ();
      fg_solve ();      
    break;

    case 3 :
      fa_setup_system ();
      fa_assemble_system ();
      fa_solve ();      
    break;

    }
}
//==================================================
//==================================================
//==================================================

void Plt_dealii::write_spherical_test(){

    if (time_step_count%time_step_solve==0) {
      switch(solve_type) {
      case 0 :
        sg_setup_system ();
        sg_assemble_system ();
        sg_solve ();      
      break;

      case 1 :
        sa_setup_system ();
        sa_assemble_system ();
        sa_solve ();      
      break;

      case 2 :
        if (re_do_setup_assemble) {
          fg_setup_system ();
          fg_assemble_system ();
        }
        fg_solve ();      
      break;

      case 3 :
        if (re_do_setup_assemble) {
          fa_setup_system ();
          fa_assemble_system ();
        }
        fa_solve ();      
      break;
      }
    }

  double c = calculate_induced_charge (0, 1);

  const auto &pos = atom_data->owned.position;      
  const auto type_i = atom_data -> owned.type [0] ;
  const auto charge_i = atom_data -> owned.charge [ type_i ];

  const dealii::Point<3> r = {pos[0].x, pos[0].y, pos[0].z};

  double p_sm = 0.0;
  dealii::Tensor<1, 3, double> f_sm_d;

  // when the point is out of the mesh, we ignore the output
  bool ignore_output = false;

  try {
    f_sm_d = - VectorTools::point_gradient (dof_handler, solution, r);
    p_sm = VectorTools::point_value (dof_handler, solution,r);
  } catch (...) {
    ignore_output = true;
  }
  
  if (ignore_output) {
  } else {
    //  Vector<double> f_sm = (f_sm_d[0], f_sm_d[1], f_sm_d[2]);
    //auto f_sm_norm = std::sqrt(f_sm*f_sm);
    auto f_sm_norm = std::sqrt(f_sm_d*f_sm_d);

    auto p_an = test_force_spherical->potential(0);
    auto f_an = test_force_spherical->field(0);
    auto f_an_norm = std::sqrt(f_an*f_an);

    auto p_si = force_field_custom[0]->potential(0);
    auto f_si = force_field_custom[0]->field(0);
    auto f_si_norm = std::sqrt(f_si*f_si);

    auto p_err = std::abs((p_si + p_sm - p_an)/p_an);
    auto f_err = std::abs((f_si_norm + f_sm_norm - f_an_norm)/f_an_norm);

    ofs_test_force_spherical << pos[0]    << " "
                             << charge_i  << " "
                             << c         << " " 
                             << p_sm      << " "
                             << f_sm_norm << " "
                             << p_si      << " "
                             << f_si_norm << " "
                             << p_an      << " "
                             << f_an_norm << " "
                             << p_si + p_sm           << " "
                             << f_si_norm + f_sm_norm << " "
                             << p_err      << " "
                             << f_err      << "\n";
;
  }

}

//==================================================
//==================================================
//==================================================


void Plt_dealii::output_vtk_solution (const std::string filename) const
{
#if defined(CAVIAR_WITH_MPI)
  if (my_mpi_rank == 0) {
#endif


  DataOut<3> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "smoothPotential");

  data_out.build_patches ();
  
  std::ofstream output (filename.c_str());
  data_out.write_vtk (output);
#if defined(CAVIAR_WITH_MPI)
  }
#endif  
}

//==================================================
//==================================================
//==================================================


void Plt_dealii::output_vtk_solution (const int cycle) const
{
#if defined(CAVIAR_WITH_MPI)
  if (my_mpi_rank == 0) {
#endif
  std::string filename ("o_solution-" +
                        Utilities::int_to_string (cycle, zeros_of_mesh_output) +
                        ".vtk");

  DataOut<3> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "smoothPotential");

  data_out.build_patches ();
  
  std::ofstream output (filename.c_str());
  data_out.write_vtk (output);
#if defined(CAVIAR_WITH_MPI)
  }
#endif  
}




//==================================================
//==================================================
//==================================================


void Plt_dealii::make_boundary_face_normals () {
  bool non_planar_face_found = false;
  face_area.clear ( );
  face_normal.clear ( );
  face_center.clear ( );
  face_id.clear ( );

  std::ofstream ofs;

#if defined(CAVIAR_WITH_MPI)
  if (my_mpi_rank == 0) 
#endif
    ofs.open("o_mesh_boundary_normals");



    for (typename Triangulation<3,3>::active_cell_iterator
         cell=triangulation.begin_active(); cell!=triangulation.end(); ++cell) {
      auto cc = cell->center();         
      
      for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f) {
        
        if (cell->face(f)->at_boundary()) {
        
          auto boundary_id = static_cast<unsigned>(cell->face(f)->boundary_id());
                              
          if (boundary_id_max < boundary_id)          
            boundary_id_max = boundary_id;  
                  
          auto fc = cell->face(f)->center();

          // n_o: a vector (approximately) to the direction of normal         
          auto n_o = cc - fc;
          
          // this is of type 'dealii::Tensor<1,3> '
          const Tensor<1,3> v01 = cell->face(f)->vertex(1) - cell->face(f)->vertex(0);
          const Tensor<1,3> v02 = cell->face(f)->vertex(2) - cell->face(f)->vertex(0);
  
          Tensor<1,3> normal = cross_product_3d(v01, v02);

          auto dot = normal * n_o;
          if (dot < 0 ) normal *= -1;
          auto norm_n = std::sqrt(normal*normal);      
          normal /= norm_n;

          //

          double area = cell->face(f)->measure ();
          if (std::isinf(area) || std::isnan(area)) {


            const Tensor<1,3> v03 = cell->face(f)->vertex(3) - cell->face(f)->vertex(0);
         

            double non_planar_measure = std::abs((v03 * normal) * (v03 * normal) /
                       ((v03 * v03) * (v01 * v01) * (v02 * v02)));
            if (non_planar_measure >=   1e-24)  {
              non_planar_face_found = true;
              //std::cout << "non planar face, measure : "<< non_planar_measure  << "\n";
            }  


            const Tensor<1,3> v12 = cell->face(f)->vertex(2) - cell->face(f)->vertex(1);          
            Tensor<1,3> twice_area = cross_product_3d(v03, v12);
            area =  0.5 * twice_area.norm();
       
          }
          //

          face_area.push_back (area);
          face_normal.push_back (normal);
          face_center.push_back (fc);
          face_id.push_back (boundary_id);
#if defined(CAVIAR_WITH_MPI)
          if (my_mpi_rank == 0)                     
#endif
            ofs << fc(0) << " " << fc(1) << " " << fc(2) << " "
                << normal[0]  << " " << normal[1]  << " " << normal[2]  << "\n";
          //std::cout << "fa " << cell->face(f)->measure () << " fi " << boundary_id << "\n";
          //std::cout << "fx:" << cell->face(f)->index() << std::endl;
        }
        
      }
    }

#if defined(CAVIAR_WITH_MPI)
  if (my_mpi_rank == 0) 
#endif
    ofs.close();
  if (non_planar_face_found)
    output->warning("Non planar face found in Plt_dealii.");

}

//==================================================
//==================================================
//==================================================

void Plt_dealii::output_boundary_id_areas () {


  std::vector<double> boundary_id_area(boundary_id_max + 1, 0);  

  for (unsigned int f = 0; f < face_id.size(); ++f) {
    boundary_id_area[face_id[f]] += face_area[f];
    //auto p1 = face_center[f];
  }

  std::ofstream ofs;

#if defined(CAVIAR_WITH_MPI)
  if (my_mpi_rank == 0) 
#endif
    ofs.open("o_mesh_boundary_area");


#if defined(CAVIAR_WITH_MPI)
  if (my_mpi_rank == 0) 
#endif
    for (unsigned int i = 0; i < boundary_id_area.size(); ++i) {
      ofs << i << " " << boundary_id_area[i] << "\n";
    }


#if defined(CAVIAR_WITH_MPI)
  if (my_mpi_rank == 0) 
#endif
    ofs.close();

}

//==================================================
//==================================================
//==================================================


double Plt_dealii::calculate_induced_charge (int t, const int requested_id) {

  std::vector<double> induced_charge(boundary_id_max + 1, 0);  



  for (unsigned int f = 0; f < face_id.size(); ++f) {

    if (face_id[f]==0) continue;
    if (std::count(face_id_ignore.begin(), face_id_ignore.end(), face_id[f]) > 0)
      continue;
    auto p1 = face_center[f];    

    //auto f_sm = - VectorTools::point_gradient (dof_handler, solution, p1);

    dealii::Tensor<1, 3, double> f_sm;
    if (ignore_point_out_of_mesh) {
      try {
        f_sm = - VectorTools::point_gradient (dof_handler, solution, p1);
      } catch (...) {
        continue;
      }
    } else {
      f_sm = - VectorTools::point_gradient (dof_handler, solution, p1);
    }

    caviar::Vector<double> pos_j {p1[0], p1[1], p1[2]};

    caviar::Vector<double> f_si {0.0, 0.0, 0.0};
    for (auto &&f_custom : force_field_custom)
      f_si += f_custom -> field(pos_j);



#ifdef CAVIAR_WITH_MPI
    double f_si_local[3] = {f_si.x, f_si.y, f_si.z};
    double f_si_total[3] = {0, 0, 0};


    MPI_Allreduce(&f_si_local, &f_si_total,
      3, MPI::DOUBLE, MPI_SUM,  MPI::COMM_WORLD);

     auto field_tot = dealii::Point<3> {f_si_total[0], f_si_total[1], f_si_total[2]} + f_sm;    
#else   
     auto field_tot = dealii::Point<3> {f_si.x, f_si.y, f_si.z} + f_sm;
#endif    
    

    auto field_normal = face_normal[f] * field_tot;

    auto local_q = field_normal * face_area[f]; // this is the charge value times 4*PI*K_el
    
    induced_charge[face_id[f]] += local_q;
    //std::cout << face_normal[f]<< " " << local_q << " " << face_area[f] << "\n";
  }

  auto four_pi_k_electrostatic_inv = 1.0/(4.0*3.14159265*k_electrostatic);

  for (auto&& i : induced_charge) {
    i *= four_pi_k_electrostatic_inv; // correcting to the absolute charge value
  }


  double sum_q = 0; // This is defined here because of CAVIAR_WITH_MPI 

#if defined(CAVIAR_WITH_MPI)
  if (my_mpi_rank == 0) {
#endif

  if (!induced_charge_id_init) {
    induced_charge_id_init = true;
    ofs_induced_charge << "# time ";
    for (unsigned int i=0; i<induced_charge.size(); ++i) {

    //if (induced_charge[i] == 0) continue; // old one
      if (std::count(face_id_ignore.begin(), face_id_ignore.end(), i) > 0)
        continue; // more general

      ofs_induced_charge << " " << i;
    }
    ofs_induced_charge << " " << "sum_q" << " " << "sum_abs_q" << "\n" << std::flush;    
  }
  
  
  // double sum_q = 0; // This is not defined here because of CAVIAR_WITH_MPI 
  double sum_abs_q = 0;  
  ofs_induced_charge << t;

  for (unsigned int i=0; i<induced_charge.size(); ++i) {
    if (std::count(face_id_ignore.begin(), face_id_ignore.end(), i) > 0)
      continue;
    auto q = induced_charge[i];
    ofs_induced_charge << " " << q;
    sum_q += q;
    sum_abs_q += std::abs(q);
  }
  ofs_induced_charge << " " << sum_q << " " << sum_abs_q << "\n" << std::flush;

#if defined(CAVIAR_WITH_MPI)
  }
#endif
  if (requested_id != -1 )
    return induced_charge[requested_id];
  return sum_q;
}


//==================================================
//==================================================
//==================================================


void Plt_dealii::run ()
{

  if (!initialized) {
    initialized = true;

    for (unsigned int i=0; i < refine_sequence_type.size(); ++i) {
      if (refine_sequence_type[i]==0)
        triangulation.refine_global (refine_sequence_value[i]);
      if (refine_sequence_type[i]==1) 
        dealii_functions::refine_boundary (triangulation,refine_sequence_value[i]);
    }
    refine_sequence_type.clear();


    make_boundary_face_normals (); 
    output_boundary_id_areas ();

    switch(solve_type) {
    case 0 :
      sg_setup_system ();
      sg_assemble_system ();
      sg_solve ();      
    break;

    case 1 :
      sa_setup_system ();
      sa_assemble_system ();
      sa_solve ();      
    break;

    case 2 :
      fg_setup_system ();
      fg_assemble_system ();
      fg_solve ();      
    break;

    case 3 :
      fa_setup_system ();
      fa_assemble_system ();
      fa_solve ();      
    break;

    }

    output_vtk_solution (0); 

    if (output_induced_charge) 
      calculate_induced_charge (0);

  } else {

    ++time_step_count;

    if (time_step_count%time_step_solve==0) {
      switch(solve_type) {
      case 0 :
        sg_setup_system ();
        sg_assemble_system ();
        sg_solve ();      
      break;

      case 1 :
        sa_setup_system ();
        sa_assemble_system ();
        sa_solve ();      
      break;

      case 2 :
        if (re_do_setup_assemble) {
          fg_setup_system ();
          fg_assemble_system ();
        }
        fg_solve ();      
      break;

      case 3 :
        if (re_do_setup_assemble) {
          fa_setup_system ();
          fa_assemble_system ();
        }
        fa_solve ();      
      break;
      }
    }     
    
    if (output_vtk && time_step_count%time_step_output_vtk==0) 
      output_vtk_solution (time_step_count);
    
    if (output_induced_charge && time_step_count%time_step_induced_charge==0) 
      calculate_induced_charge (time_step_count);
  }

}

//==================================================
//==================================================
//==================================================


void Plt_dealii::run_time_profile ()
{

  if (!initialized) {
    initialized = true;

    for (unsigned int i=0; i < refine_sequence_type.size(); ++i) {
      if (refine_sequence_type[i]==0)
        triangulation.refine_global (refine_sequence_value[i]);
      if (refine_sequence_type[i]==1) 
        dealii_functions::refine_boundary (triangulation,refine_sequence_value[i]);
    }
    refine_sequence_type.clear();


    make_boundary_face_normals (); 
    output_boundary_id_areas ();

    clock_t t1=0, t2=0, t3=0, t4=0, t5=0, t6=0;
    switch(solve_type) {
    case 0 :
      t1 = clock();
      sg_setup_system ();
      t2 = clock();
      sg_assemble_system ();
      t3 = clock();
      sg_solve ();      
      t4 = clock();
    break;

    case 1 :
      t1 = clock();
      sa_setup_system ();
      t2 = clock();
      sa_assemble_system ();
      t3 = clock();
      sa_solve ();      
      t4 = clock();
    break;

    case 2 :
      t1 = clock();
      fg_setup_system ();
      t2 = clock();
      fg_assemble_system ();
      t3 = clock();
      fg_solve ();      
      t4 = clock();
    break;

    case 3 :
      t1 = clock();
      fa_setup_system ();
      t2 = clock();
      fa_assemble_system ();
      t3 = clock();
      fa_solve_time_profile ();      
      t4 = clock();
    break;

    }

    output_vtk_solution (0); 

    t5 = clock();
    if (output_induced_charge) 
      calculate_induced_charge (0);
    t6 = clock();

    std::cout << "setup: "    << (double)(t2 - t1)/CLOCKS_PER_SEC << "\n";
    std::cout << "assemble: " << (double)(t3 - t2)/CLOCKS_PER_SEC << "\n";
    std::cout << "solve: "    << (double)(t4 - t3)/CLOCKS_PER_SEC << "\n";
    std::cout << "induced_charge: " << (double)(t6 - t5)/CLOCKS_PER_SEC << "\n";


  } else {

    ++time_step_count;
      clock_t t1=0, t2=0, t3=0, t4=0, t5=0, t6=0;
    if (time_step_count%time_step_solve==0) {

      t1 = clock();
      t2 = clock();
      switch(solve_type) {
      case 0 :
        t1 = clock();
        sg_setup_system ();
        t2 = clock();
        sg_assemble_system ();
        t3 = clock();
        sg_solve ();      
        t4 = clock();
      break;

      case 1 :
        t1 = clock();
        sa_setup_system ();
        t2 = clock();
        sa_assemble_system ();
        t3 = clock();
        sa_solve ();      
        t4 = clock();
      break;

      case 2 :
        if (re_do_setup_assemble) {
          t1 = clock();
          fg_setup_system ();
          t2 = clock();
          fg_assemble_system ();
          t3 = clock();
        }
        t3 = clock();
        fg_solve ();      
        t4 = clock();
      break;

      case 3 :
        if (re_do_setup_assemble) {
          t1 = clock();
          fa_setup_system ();
          t2 = clock();
          fa_assemble_system ();
          t3 = clock();
        }
        t3 = clock();
        fa_solve_time_profile ();      
        t4 = clock();
      break;
      }
    }     
    
    if (output_vtk && time_step_count%time_step_output_vtk==0) 
      output_vtk_solution (time_step_count);

    t5 = clock();
    if (output_induced_charge && time_step_count%time_step_induced_charge==0) 
      calculate_induced_charge (time_step_count);
    t6 = clock();

    std::cout << "setup: "    << (double)(t2 - t1)/CLOCKS_PER_SEC << "\n";
    std::cout << "assemble: " << (double)(t3 - t2)/CLOCKS_PER_SEC << "\n";
    std::cout << "solve: "    << (double)(t4 - t3)/CLOCKS_PER_SEC << "\n";
    std::cout << "induced_charge: " << (double)(t6 - t5)/CLOCKS_PER_SEC << "\n";
  }

}

//==================================================
//==================================================
//==================================================

void Plt_dealii::verify_settings() {
  FC_NULLPTR_CHECK(atom_data)
}

//==================================================
//==================================================
//==================================================


void Plt_dealii::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS

  if (time_profile) {
    clock_t t1 = clock();

    run_time_profile ();
    clock_t t2 = clock();

    calculate_all_particles_mesh_force_acc();

    clock_t t3 = clock();

    double t2_1= (double)(t2 - t1)/CLOCKS_PER_SEC;
    double t3_2= (double)(t3 - t2)/CLOCKS_PER_SEC;

    std::cout << "total_run      : " << t2_1 << "\n";
    std::cout << "mesh_force_acc : " << t3_2 << "\n";
  } else {
    run ();

    calculate_all_particles_mesh_force_acc();
  }

#ifdef CAVIAR_WITH_MPI
  // XXX: do not delete this part yet!
  // In case of 4 threads, this thing happens
  //
  // 0 -> 1, 2, 3, 4 -> 0
  // 1 -> 2, 3, 4, -> 1
  // 2 -> 3, 4 -> 2
  // 3 -> 4 -> 3
  /*
  error->all(FC_FILE_LINE_FUNC,"parallel code not implemented yet");

  int me, nprocs;
  MPI_Comm_rank (MPI::COMM_WORLD, &me);
  MPI_Comm_size (MPI::COMM_WORLD, &nprocs);

  MPI_Barrier (MPI_COMM_WORLD);  
  const auto &pos = atom_data -> owned.position;
  const auto &type = atom_data -> owned.type;    
  int root = 0;
  while (root < nprocs) {
  
    unsigned root_pos_size = pos.size();
    MPI_Bcast (&root_pos_size, 1,  MPI::UNSIGNED, root, MPI::COMM_WORLD);
    MPI_Barrier (MPI::COMM_WORLD);
      
    if (root == me) {
      

      for (int i = root+1; i < nprocs; ++i) {
        MPI_Send (type.data(), root_pos_size, MPI::UNSIGNED, i, 0, MPI::COMM_WORLD);
        MPI_Send (pos.data(), 3*root_pos_size, MPI::DOUBLE, i, 1, MPI::COMM_WORLD);  
      }

      std::vector<Vector<double>> root_acc (root_pos_size,{0,0,0});
      
      for (int i = root+1; i < nprocs; ++i) {
        MPI_Recv (root_acc.data(), 3*root_pos_size, MPI::DOUBLE,
          root, 2, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
      }
      
      for (unsigned  i=0;i<pos.size();++i) {      
        atom_data -> owned.acceleration [i] += root_acc[i];
      }
      
    } else {
      if (me > root) {
        std::vector<Vector<double>> root_pos   (root_pos_size);
        std::vector<unsigned> root_type (root_pos_size);        
        
        MPI_Recv (root_type.data(), root_pos_size, MPI::UNSIGNED,
          root, 0, MPI::COMM_WORLD, MPI_STATUS_IGNORE);        
          
        MPI_Recv (root_pos.data(), 3*root_pos_size, MPI::DOUBLE,
          root, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);        


        std::vector<Vector<double>> root_acc (root_pos_size,{0,0,0});    


        for (unsigned i=0;i<pos.size();++i) {
           const auto type_i = atom_data -> owned.type [i] ;
          const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];
          const auto charge_i = atom_data -> owned.charge [ type_i ];

   
          for (unsigned j=i+1;j<root_pos.size();++j) {    
             const auto type_j = root_type [j] ;
            const auto mass_inv_j = atom_data -> owned.mass_inv [ type_j ];
            const auto charge_j = atom_data -> owned.charge [ type_j ];      
            const auto dr = root_pos[j] - pos[i]; 
            const auto dr_sq = dr*dr;
            const auto dr_norm = std::sqrt(dr_sq);      
            const auto force = k_electrostatic * charge_i * charge_j * dr / (dr_sq*dr_norm);
            atom_data -> owned.acceleration [i] -= force * mass_inv_i;
            root_acc [j] += force * mass_inv_j; // no periodic boundary condition yet      
          }
        }
      }      
    }
    
    MPI_Barrier (MPI::COMM_WORLD);  
    ++root;
  }  
  */  
  
#endif


}

//==================================================
//==================================================
//==================================================

void Plt_dealii::calculate_all_particles_mesh_force_acc() {
  const auto &pos = atom_data -> owned.position;    
  for (unsigned int i=0;i<pos.size();++i) {
    const auto type_i = atom_data -> owned.type [i] ;
    const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];
    const auto charge_i = atom_data -> owned.charge [ type_i ];

    const dealii::Point<3> r = {pos[i].x, pos[i].y, pos[i].z};
    dealii::Tensor<1, 3, double> field;
    if (ignore_point_out_of_mesh) {
      try {
        field = - VectorTools::point_gradient (dof_handler, solution, r);
      } catch (...) {
        continue;
      }
    } else {
      field = - VectorTools::point_gradient (dof_handler, solution, r);
    }
    auto frc = field * charge_i;
    auto force = caviar::Vector<double> {frc[0], frc[1], frc[2]};
    atom_data -> owned.acceleration [i] += force * mass_inv_i;    
  }
}

//==================================================
//==================================================
//==================================================


} //force_field
} //objects
} // namespace caviar

#else

namespace caviar {
namespace objects {
namespace force_field {

Plt_dealii::Plt_dealii(CAVIAR * fptr) : 
    caviar::objects::Force_field{fptr} { 
  error->all(FC_FILE_LINE_FUNC,"please recompile with DEAL.II library.");
}

Plt_dealii::~Plt_dealii() {}

bool Plt_dealii::read (caviar::interpreter::Parser *parser) {
  error->all(FC_FILE_LINE_FUNC,"please recompile with DEAL.II library.");
  parser->get_val_token();
  return true;
}

void Plt_dealii::calculate_acceleration () {
  error->all(FC_FILE_LINE_FUNC,"please recompile with DEAL.II library.");
}

} //force_field
} //objects
} // namespace caviar
#endif
