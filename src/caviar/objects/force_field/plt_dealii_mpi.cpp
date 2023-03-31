
//========================================================================
//
// Copyright (C) 2019 by deal.II authors and
// Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Some part of this file (Solving Laplace equation) is based on the step-40
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

#include "caviar/objects/force_field/plt_dealii_mpi.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/interpreter/communicator.h"

#ifdef CAVIAR_WITH_DEALII_MPI


#include "caviar/objects/atom_data.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/utility/macro_constants.h"

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


CAVIAR_NAMESPACE_OPEN

namespace force_field {

//==================================================
//==================================================
//==================================================



Plt_dealii_mpi::Plt_dealii_mpi(CAVIAR * fptr) : 
    caviar::Force_field{fptr},
    triangulation (mpi_comm,
                   typename Triangulation<3>::MeshSmoothing
                   (Triangulation<3>::smoothing_on_refinement |
                    Triangulation<3>::smoothing_on_coarsening)),
    tria_reserve (mpi_comm,
                   typename Triangulation<3>::MeshSmoothing
                   (Triangulation<3>::smoothing_on_refinement |
                    Triangulation<3>::smoothing_on_coarsening)),
    fe (FC_DEALII_FE_Q_POLY_DEGREE),
    dof_handler (triangulation),
    num_quadrature_points{2}, 
    solver_control_tolerance{1e-6},
    use_preconditioner{true},
    atom_data{nullptr}, 
    initialized{false}, k_electrostatic{1.0},
    time_step_count{0}, time_step_solve{1},
    time_step_output_vtk{1}, time_step_induced_charge{1},
    output_vtk{false},
    induced_charge_id_init{false},
    output_induced_charge{false}, boundary_id_max{0},
    make_time_profile{false},
    zeros_of_mesh_output{2},
    my_mpi_rank{fptr->comm->me},
    mpi_world_size{fptr->comm->nprocs}
{ 
  FC_OBJECT_INITIALIZE_INFO

  ignore_point_out_of_mesh = false;
  output->warning( " "
      "Always add force_field::DealII_poisson_custom_mpi after 'ewald_k' and"
      " 'slab' to the integrators because they have to be initialized in"
      " every steps. The initialization functions are done in "
      "calculate_acceleration function.");


}

//==================================================
//==================================================
//==================================================

Plt_dealii_mpi::~Plt_dealii_mpi() {
    dof_handler.clear ();
}

//==================================================
//==================================================
//==================================================


bool Plt_dealii_mpi::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"refine_global")) {    
      int n = parser->get_literal_int();     
      refine_sequence_type.push_back (0);
      refine_sequence_value.push_back (n);
    } else if (string_cmp(t,"num_quadrature_points")) {    
      num_quadrature_points = parser->get_literal_int();     
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
    } else if (string_cmp(t,"solver_control_tolerance")) {
      solver_control_tolerance = parser->get_literal_real();
    } else if (string_cmp(t,"add_unv_mesh")) {
      auto token = parser->get_val_token();
      auto file_name = token.string_value;
      unv_mesh_filename.push_back(file_name);
    } else if (string_cmp(t,"read_unv_mesh")) {
      //error->all (FC_FILE_LINE_FUNC_PARSE, "not implemented yet");
      dealii_functions::read_domain(triangulation, tria_reserve, unv_mesh_filename);
    } else if (string_cmp(t,"time_step_solve")) {
      time_step_solve = parser->get_literal_int();    
    } else if (string_cmp(t,"output_vtk")) {
      output_vtk = true;
      time_step_output_vtk = parser->get_literal_int();
    } else if (string_cmp(t,"output_induced_charge")) {
      output_induced_charge = true;
      time_step_induced_charge = parser->get_literal_int();
      if (my_mpi_rank == 0)     
        ofs_induced_charge.open ("o_induced_charge");

    } else if (string_cmp(t,"induced_charge_ignore_id")) {
      unsigned id = parser->get_literal_int();
      face_id_ignore.push_back(id);
    } else if (string_cmp(t,"dealii_grid_generator")) {
      dealii_functions::dealii_grid_generator(fptr, parser, triangulation);
      return in_file;
    } else if (string_cmp(t,"add_force_field") || string_cmp(t,"force_field")) {
      FIND_OBJECT_BY_NAME(force_field,it)
      force_field_custom.push_back(object_container->force_field[it->second.index]);
    } else if (string_cmp(t,"test_code")) {
      test_code();
    } else if (string_cmp(t,"ignore_point_out_of_mesh")) {
      ignore_point_out_of_mesh = true;
    } else if (string_cmp(t,"re_initialize")) {
      initialized = false;
    } else if (string_cmp(t,"calculate_acceleration")) {
      calculate_acceleration();
    }
    
    else FC_ERR_UNDEFINED_VAR(t)
  }
  
  return in_file;
}

//==================================================
//==================================================
//==================================================
void Plt_dealii_mpi::test_code() {

/*
  atom_data->synch_owned_data(0);

  //GridGenerator::hyper_cube (triangulation);
  //triangulation.refine_global (5);
  std::string st = "../data/mesh_periodic_fixed.unv";
  std::ifstream in;
  in.open(st.c_str());
  GridIn<3,3> gi;
  gi.attach_triangulation(triangulation);
  gi.read_unv(in);


  //dealii_functions::read_domain(triangulation, tria_reserve, unv_mesh_filename);

  setup_system ();

  assemble_system ();

  solve ();

  output_vtk_solution (0); 

  calculate_all_particles_mesh_force_acc();

  calculate_induced_charge (0);*/

}

//==================================================
//==================================================
//==================================================


void Plt_dealii_mpi::output_vtk_solution (const int cycle) const
{
    DataOut<3> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (locally_relevant_solution, "smoothPotential");

    dealii::Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");

    data_out.build_patches ();

    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (cycle, zeros_of_mesh_output) +
                                  "." +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);

    if (Utilities::MPI::this_mpi_process(mpi_comm) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i=0;
             i<Utilities::MPI::n_mpi_processes(mpi_comm);
             ++i)
          filenames.push_back ("solution-" +
                               Utilities::int_to_string (cycle, zeros_of_mesh_output) +
                               "." +
                               Utilities::int_to_string (i, 4) +
                               ".vtu");

        std::ofstream master_output (("solution-" +
                                      Utilities::int_to_string (cycle, 2) +
                                      ".pvtu").c_str());
        data_out.write_pvtu_record (master_output, filenames);
      }
}

//==================================================
//==================================================
//==================================================


void Plt_dealii_mpi::make_boundary_face_normals () {
  bool non_planar_face_found = false;
  face_area.clear ( );
  face_normal.clear ( );
  face_center.clear ( );
  face_id.clear ( );

  std::string fn = "o_mesh_boundary_normals-" + std::to_string(my_mpi_rank);
  std::ofstream ofs (fn.c_str()); 


  // for (typename Triangulation<3,3>::active_cell_iterator
  //   cell=triangulation.begin_active(); cell!=triangulation.end(); ++cell)
  typename DoFHandler<3>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())  {

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
                    
          ofs << fc(0) << " " << fc(1) << " " << fc(2) << " "
              << normal[0]  << " " << normal[1]  << " " << normal[2]  << "\n";
              

        }
        
      }
    }

  if (non_planar_face_found)
    output->warning("Non planar face found in Plt_dealii.");

}

//==================================================
//==================================================
//==================================================


void Plt_dealii_mpi::calculate_induced_charge (int t) {


  MPI_Allreduce (MPI::IN_PLACE, &boundary_id_max, 1, MPI::INT, MPI_MAX, mpi_comm);

  std::vector<double> induced_charge(boundary_id_max + 1, 0);  

  for (unsigned int f = 0; f < face_id.size(); ++f) {

    if (face_id[f]==0) continue;

    if (std::count(face_id_ignore.begin(), face_id_ignore.end(), face_id[f]) > 0)
      continue;

    auto p1 = face_center[f];    

    
    // smooth field
    dealii::Tensor<1, 3, double> f_sm;
    f_sm = - VectorTools::point_gradient (dof_handler, locally_relevant_solution, p1); 


#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

    caviar::Vector<double> pos_j {p1[0], p1[1], p1[2]};

    // singular field
    caviar::Vector<double> f_si {0.0, 0.0, 0.0};

    for (auto &&f_custom : force_field_custom)
      f_si += f_custom -> field(pos_j);

    
    auto field_tot = dealii::Point<3> {f_si.x, f_si.y, f_si.z} + f_sm;          

    auto field_normal = face_normal[f] * field_tot;

    auto local_q = field_normal * face_area[f];  // this is the charge value times 4*PI*K_el

    induced_charge[face_id[f]] += local_q;

#else
    error -> all (FC_FILE_LINE_FUNC, "This part is not implemented.");
#endif
  }

  auto four_pi_k_electrostatic_inv = 1.0/(4.0*3.14159265*k_electrostatic);

  for (auto&& i : induced_charge) {
    i *= four_pi_k_electrostatic_inv; // correcting to the absolute charge value
  }

  // used for debugging. Output charges by process.
  /*
  static std::string ofs_st = "ind-" + std::to_string(my_mpi_rank);
  static std::ofstream ofs(ofs_st.c_str());

  ofs << t ;
  for (auto&& q : induced_charge) {
      ofs << " " << q;
  }
  ofs << std::endl;
  */

  if (mpi_world_size!=1) 
    MPI_Allreduce (MPI::IN_PLACE, &induced_charge[0], boundary_id_max + 1, MPI::DOUBLE, MPI_SUM, mpi_comm);


  if (my_mpi_rank == 0) {

    if (!induced_charge_id_init) {
      induced_charge_id_init = true;
      ofs_induced_charge << "# time ";
      for (unsigned int i=0; i<induced_charge.size(); ++i) {
        if (induced_charge[i] == 0) continue;
        ofs_induced_charge << " " << i;
      }
    ofs_induced_charge << " " << "sum_q" << " " << "sum_abs_q" << "\n" << std::flush;    
    }
  
    double sum_q = 0;
    double sum_abs_q = 0;  
    ofs_induced_charge << t;
    for (auto&& q : induced_charge) {
      if (q == 0) continue;
      ofs_induced_charge << " " << q;
      sum_q += q;
      sum_abs_q += std::abs(q);
    }
    ofs_induced_charge << " " << sum_q << " " << sum_abs_q << "\n" << std::flush;
  }


}
//==================================================
//==================================================
//==================================================


void Plt_dealii_mpi::run ()
{
  if (!initialized) {

    std::string st_com = "Running with ";
#ifdef USE_PETSC_LA
    st_com += "PETSc";
#else
    st_com += "Trilinos";
#endif
    st_com += " on " + std::to_string( Utilities::MPI::n_mpi_processes(mpi_comm));
    st_com += " MPI rank(s).";
    output->info(st_com);


    initialized = true;


    for (unsigned int i=0; i < refine_sequence_type.size(); ++i) {
      output->warning( " "
        "Refining a distributed mesh is may affect boundary ids."
        " ");

      if (refine_sequence_type[i]==0)
        triangulation.refine_global (refine_sequence_value[i]);
      if (refine_sequence_type[i]==1) 
        dealii_functions::refine_boundary (triangulation, refine_sequence_value[i]);
    }
    refine_sequence_type.clear();


    make_boundary_face_normals ();


    sa_setup_system ();

    sa_assemble_system ();

    sa_solve ();      


    output_vtk_solution (0); 

    if (output_induced_charge) 
      calculate_induced_charge (0);

  } else {


    ++time_step_count;


    if (time_step_count%time_step_solve==0) {
      sa_setup_system ();
      sa_assemble_system ();
      sa_solve ();
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

void Plt_dealii_mpi::verify_settings() {
#ifndef CAVIAR_SINGLE_MPI_MD_DOMAIN
  error->all(FC_FILE_LINE_FUNC,"This configuration is not implemented for this "
      "force_field yet, meaning, it cannot handle distributed mesh alongside "
      "distributed MD domains. You may re-compile CAVIAR with "
      "CAVIAR_SINGLE_MPI_MD_DOMAIN option, or use Plt_dealii with "
      "with mpirun, meaning, all distributed particles share the same mesh." );
#endif

  FC_NULLPTR_CHECK(atom_data)
}

//==================================================
//==================================================
//==================================================


void Plt_dealii_mpi::calculate_acceleration () {
  FC_OBJECT_VERIFY_SETTINGS


  run ();
  
  calculate_all_particles_mesh_force_acc();
  

  
}


//==================================================
//==================================================
//==================================================

void Plt_dealii_mpi::calculate_all_particles_mesh_force_acc() {

 
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

  /* // DATA check
  auto &pos = atom_data->owned.position;  
  caviar::Vector<double> pm{0,0,0};
  for (unsigned int i=0;i<pos.size();++i) {
    pm += pos[i];
  }
  pm /= pos.size();
  std::cout << my_mpi_rank << " , pm : " << pm << std::endl;
  */

  auto mpi_fc_vector_type = comm->mpi_fc_vector_type;
  auto &pos = atom_data->owned.position;
  auto pos_size = pos.size();

  static std::vector<int> pos_in_mesh(pos_size, 0);

  static std::vector<caviar::Vector<Real_t>> acc_found(pos_size, caviar::Vector<Real_t> {0,0,0});

  

  // additional acceleration calculation for each domain

  if (my_mpi_rank == 0) {

    for (unsigned int i=0;i<pos.size();++i) {

      const dealii::Point<3> r = {pos[i].x, pos[i].y, pos[i].z};
      dealii::Tensor<1, 3, double> field;

      try {
        field = - VectorTools::point_gradient (dof_handler, locally_relevant_solution, r);

      } catch (...) {
        pos_in_mesh[i]=false;
        continue;
      }
      pos_in_mesh[i]=true;

      const auto type_i = atom_data -> owned.type [i] ;
      const auto mass_inv_i = atom_data -> owned.mass_inv [ type_i ];
      const auto charge_i = atom_data -> owned.charge [ type_i ];

      auto frc = field * charge_i;
      auto force = caviar::Vector<double> {frc[0], frc[1], frc[2]};

      atom_data -> owned.acceleration [i] += force * mass_inv_i;

    }

  } else {

    for (unsigned int i=0;i<pos.size();++i) {

      const dealii::Point<3> r = {pos[i].x, pos[i].y, pos[i].z};
      dealii::Tensor<1, 3, double> field;

      try {
        field = - VectorTools::point_gradient (dof_handler, locally_relevant_solution, r);
      } 
      catch (std::exception &exc) {
        pos_in_mesh[i] = false;
        acc_found[i] = caviar::Vector<double> {0.,0.,0.};
        continue;
      }
      catch (...) {
        pos_in_mesh[i] = false;
        acc_found[i] = caviar::Vector<double> {0.,0.,0.};
        continue;
      }
      pos_in_mesh[i]=true;

      const auto type_i   = atom_data -> owned.type [i] ;
      const auto mass_inv_i   = atom_data -> owned.mass_inv [ type_i ];
      const auto charge_i = atom_data -> owned.charge [ type_i ];

      auto frc = field * charge_i;
      auto force = caviar::Vector<double> {frc[0], frc[1], frc[2]};

      acc_found[i] = force* mass_inv_i;

    }
  }


  if (mpi_world_size==1) return;

  // sending additional acceleration to the root process and adding it.
  if (my_mpi_rank == 0) {

    for (int i = 1; i < mpi_world_size; ++i) {

      MPI_Recv (&pos_in_mesh[0],   pos_size, MPI_INT,            i, 0, mpi_comm, MPI_STATUS_IGNORE); 
      MPI_Recv (&acc_found[0],     pos_size, mpi_fc_vector_type, i, 1, mpi_comm, MPI_STATUS_IGNORE); 

      for (unsigned int j = 0; j < pos_size; ++j) {
        if (pos_in_mesh[j]==true) {
          atom_data -> owned.acceleration [j] += acc_found[j];
        }
      }


    }


  } else {
    MPI_Send (&pos_in_mesh[0],   pos_size, MPI_INT,            0, 0, mpi_comm); 
    MPI_Send (&acc_found[0],     pos_size, mpi_fc_vector_type, 0, 1, mpi_comm); 
  }


#else

#endif

}

//==================================================
//==================================================
//==================================================


} //force_field

CAVIAR_NAMESPACE_CLOSE

#else

CAVIAR_NAMESPACE_OPEN

namespace force_field {

Plt_dealii_mpi::Plt_dealii_mpi(CAVIAR * fptr) : 
    caviar::Force_field{fptr} {
  error->all(FC_FILE_LINE_FUNC,"please recompile with DEAL.II library.");
}

Plt_dealii_mpi::~Plt_dealii_mpi() {}

bool Plt_dealii_mpi::read (caviar::interpreter::Parser *parser) {
  error->all(FC_FILE_LINE_FUNC,"please recompile with DEAL.II library.");
  parser->get_val_token();
  return true;
}

void Plt_dealii_mpi::calculate_acceleration () {
  error->all(FC_FILE_LINE_FUNC,"please recompile with DEAL.II library.");
}

} //force_field

CAVIAR_NAMESPACE_CLOSE
#endif
