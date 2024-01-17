
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
#include "caviar/objects/domain.h"
#include "caviar/utility/macro_constants.h"
#include "caviar/utility/time_utility.h"
#include "caviar/utility/file_utility.h"
#include "caviar/objects/unique/time_function.h"
#include "caviar/objects/unique/time_function_3d.h"
#include "caviar/objects/unique/grid_1d.h"

#ifdef CAVIAR_WITH_MPI
#include <mpi.h>
#endif

#ifdef CAVIAR_WITH_DEALII

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

// #include <deal.II/grid/grid_all_interpreter_tools.h>
#include <deal.II/grid/grid_tools.h>

#if DEALII_VERSION_MAJOR == 8
#include <deal.II/grid/tria_boundary_lib.h>
#endif

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

// #include <deal.II/dofs/dof_all_interpreter_tools.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

// #include <deal.II/numerics/vector_all_interpreter_tools.h>
#include <deal.II/numerics/vector_tools.h>

// #include <deal.II/numerics/matrix_all_interpreter_tools.h>
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
#if DEALII_VERSION_MAJOR == 8
#include <deal.II/lac/filtered_matrix.h>
#endif

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
#include "caviar/objects/neighborlist.h" // used for ml training

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  //==================================================
  //==================================================
  //==================================================

  Plt_dealii::Plt_dealii(CAVIAR *fptr) : caviar::Force_field{fptr},
                                         fe(FC_DEALII_FE_Q_POLY_DEGREE),
                                         dof_handler(triangulation),
                                         num_quadrature_points{FC_DEALII_FE_Q_POLY_DEGREE + 1},
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
                                         ,
                                         my_mpi_rank{fptr->comm->me}, mpi_world_size{fptr->comm->nprocs}
#endif
  {
    FC_OBJECT_INITIALIZE_INFO
    ignore_point_out_of_mesh = false;
    time_profile = false;
    output->warning(" "
                    "Always add force_field::DealII_poisson_custom after 'ewald_k' and"
                    " 'slab' to the integrators because they have to be initialized in"
                    " every steps. The initialization functions are done in "
                    "calculate_acceleration function.");
    test_force_spherical = nullptr;
    init_test_force_spherical = false;
  }

  //==================================================
  //==================================================
  //==================================================

  Plt_dealii::~Plt_dealii() {}

  //==================================================
  //==================================================
  //==================================================

  bool Plt_dealii::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR
      if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else if (string_cmp(t, "make_grid"))
      {
        // make_grid();
        error->all(FC_FILE_LINE_FUNC_PARSE, "not implemented");
      }
      else if (string_cmp(t, "num_quadrature_points"))
      {
        num_quadrature_points = parser->get_literal_int();
      }
      else if (string_cmp(t, "refine_global"))
      {
        int n = parser->get_literal_int();
        refine_sequence_type.push_back(0);
        refine_sequence_value.push_back(n);
      }
      else if (string_cmp(t, "refine_boundary"))
      {
        int n = parser->get_literal_int();
        refine_sequence_type.push_back(1);
        refine_sequence_value.push_back(n);
      }
      else if (string_cmp(t, "boundary_id_value"))
      {
        int id = parser->get_literal_int();
        if (id == 0)
        {
          error->all(FC_FILE_LINE_FUNC_PARSE, "#boundary_id=0 is reserved for free boundary");
        }
        double value = 0;
        GET_OR_CHOOSE_A_REAL(value, "", "")
        boundary_id_value.push_back(std::make_pair(id, value));
      }
      else if (string_cmp(t, "boundary_id_time_function"))
      {
        int id = parser->get_literal_int();
        if (id == 0)
        {
          error->all(FC_FILE_LINE_FUNC_PARSE, "#boundary_id=0 is reserved for free boundary");
        }
        FIND_OBJECT_BY_NAME(unique, it)
        FC_CHECK_OBJECT_CLASS_NAME(unique, it, time_function)
        unique::Time_function *a = dynamic_cast<unique::Time_function *>(object_container->unique[it->second.index]);
        boundary_id_time_function.push_back(std::make_pair(id, a));
      }
      else if (string_cmp(t, "k_electrostatic"))
      {
        GET_OR_CHOOSE_A_REAL(k_electrostatic, "", "")
        if (k_electrostatic < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "k_electrostatic has to be non-negative.");
      }
      else if (string_cmp(t, "preconditioner_relaxation"))
      {
        preconditioner_relaxation = parser->get_literal_real();
        if (preconditioner_relaxation < 0 || preconditioner_relaxation > 2)
          error->all(FC_FILE_LINE_FUNC, "expected '0.0 < preconditioner_relaxation < 2.0'.");
      }
      else if (string_cmp(t, "solver_control_tolerance"))
      {
        solver_control_tolerance = parser->get_literal_real();
      }
      else if (string_cmp(t, "add_unv_mesh"))
      {
        auto token = parser->get_val_token();
        auto file_name = token.string_value;
        std::string file_name_full = join_path(fptr->input_file_directory, file_name);
        if (!file_exists_1(file_name_full))
          error->all(FC_FILE_LINE_FUNC_PARSE, "file does not exist : " + file_name_full);
        unv_mesh_filename.push_back(file_name_full);
      }
      else if (string_cmp(t, "read_unv_mesh"))
      {
        dealii_functions::read_domain(triangulation, tria_reserve, unv_mesh_filename);
      }
      else if (string_cmp(t, "set_domain") || string_cmp(t, "domain"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
      }
      else if (string_cmp(t, "time_step_solve"))
      {
        time_step_solve = parser->get_literal_int();
      }
      else if (string_cmp(t, "output_vtk"))
      {
        output_vtk = true;
        time_step_output_vtk = parser->get_literal_int();
      }
      else if (string_cmp(t, "output_induced_charge"))
      {
        output_induced_charge = true;
        time_step_induced_charge = parser->get_literal_int();
#if defined(CAVIAR_WITH_MPI)
        if (my_mpi_rank == 0)
          ofs_induced_charge.open("o_induced_charge");
#else
        ofs_induced_charge.open("o_induced_charge");
#endif
      }
      else if (string_cmp(t, "output_induced_charge_name"))
      {
        output_induced_charge = true;
        time_step_induced_charge = parser->get_literal_int();
        auto token = parser->get_val_token();
        std::string fn = token.string_value;
#if defined(CAVIAR_WITH_MPI)
        if (my_mpi_rank == 0)
          ofs_induced_charge.open(fn.c_str());
#else
        ofs_induced_charge.open(fn.c_str());
#endif
      }
      else if (string_cmp(t, "induced_charge_ignore_id"))
      {
        unsigned id = parser->get_literal_int();
        face_id_ignore.push_back(id);
      }
      else if (string_cmp(t, "dealii_grid_generator"))
      {
        dealii_functions::dealii_grid_generator(fptr, parser, triangulation);
        return in_file;
      }
      else if (string_cmp(t, "add_force_field") || string_cmp(t, "force_field"))
      {
        FIND_OBJECT_BY_NAME(force_field, it)
        force_field_custom.push_back(object_container->force_field[it->second.index]);
      }
      else if (string_cmp(t, "set_solve_type"))
      {
        auto t1 = parser->get_val_token();
        auto st = t1.string_value;
        if (string_cmp(st, "simple_global"))
        {
          solve_type = 0;
        }
        else if (string_cmp(st, "simple_adaptive"))
        {
          solve_type = 1;
        }
        else if (string_cmp(st, "faster_global"))
        {
          solve_type = 2;
        }
        else if (string_cmp(st, "faster_adaptive"))
        {
          solve_type = 3;
        }
        else
        {
          error->all(FC_FILE_LINE_FUNC, static_cast<std::string>("undefined option : ") + st);
        }
      }
      else if (string_cmp(t, "re_initialize"))
      {
        initialized = false;
      }
      else if (string_cmp(t, "calculate_acceleration"))
      {
        calculate_acceleration();
      }
      else if (string_cmp(t, "ignore_point_out_of_mesh"))
      {
        ignore_point_out_of_mesh = true;
      }
      else if (string_cmp(t, "time_profile"))
      {
        time_profile = true;
      }
      else if (string_cmp(t, "write_spherical_test"))
      {
        write_spherical_test();
      }
      else if (string_cmp(t, "set_spherical_test_force"))
      {
        FIND_OBJECT_BY_NAME(force_field, it)
        test_force_spherical = object_container->force_field[it->second.index];
      }
      else if (string_cmp(t, "start_spherical_test"))
      {
        auto token = parser->get_val_token();
        spherical_test_file_name = token.string_value;
        ofs_test_force_spherical.open(spherical_test_file_name.c_str());
        start_spherical_test();
      }
      else if (string_cmp(t, "finish_spherical_test"))
      {
        ofs_test_force_spherical.close();
        ofs_induced_charge.close();
      }
      else if (string_cmp(t, "refine_global_here"))
      {
        triangulation.refine_global(1);
      }
      else if (string_cmp(t, "output_vtk_here"))
      {
        auto token = parser->get_val_token();
        spherical_test_file_name = token.string_value;
        output_vtk_solution(spherical_test_file_name);
      }
      else if (string_cmp(t, "output_field_vectors"))
      {
        output_field_vectors(parser);
        return in_file;
      }
      else if (string_cmp(t, "output_potential_values"))
      {
        output_potential_values(parser);
        return in_file;
      }
      else if (string_cmp(t, "generate_ml_training_data"))
      {
        generate_ml_training_data(parser);
        return in_file;
      }
      else if (string_cmp(t, "set_position_offset"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        FC_CHECK_OBJECT_CLASS_NAME(unique, it, time_function_3d)
        unique::Time_function_3d *a = dynamic_cast<unique::Time_function_3d *>(object_container->unique[it->second.index]);
        position_offset = a;
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::generate_ml_training_data(caviar::interpreter::Parser *parser)
  {

    unique::Grid_1D *grid_1d_x = nullptr, *grid_1d_y = nullptr, *grid_1d_z = nullptr;
    Neighborlist *neighborlist = nullptr;
    bool in_file = true;
    if (in_file == true)
    {
      // removes a warning
    }
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION

      auto t = token.string_value;

      if (string_cmp(t, "neighborlist"))
      {
        FIND_OBJECT_BY_NAME(neighborlist, it)
        neighborlist = static_cast<Neighborlist *>(object_container->neighborlist[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_x"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_1d_x = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_y"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_1d_y = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_z"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_1d_z = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }

      else
        FC_ERR_UNDEFINED_VAR(t)
    }
    FC_NULLPTR_CHECK(neighborlist)
    FC_NULLPTR_CHECK(grid_1d_x)
    FC_NULLPTR_CHECK(grid_1d_y)
    FC_NULLPTR_CHECK(grid_1d_z)
    FC_NULLPTR_CHECK(atom_data)

    std::cout << "generate_ml_training_data: Check condition " << std::endl;
    // check condition
    auto &pos = atom_data->atom_struct_owned.position;
    if (pos.size() != 1)
      error->all(FC_FILE_LINE_FUNC, "Expected only one atom in atomdata for ML training data.");

    std::cout << "generate_ml_training_data: init FE calculation" << std::endl;
    // init FE calculation
    if (!initialized)
    {
      initialized = true;

      for (unsigned int i = 0; i < refine_sequence_type.size(); ++i)
      {
        if (refine_sequence_type[i] == 0)
          triangulation.refine_global(refine_sequence_value[i]);
        if (refine_sequence_type[i] == 1)
          dealii_functions::refine_boundary(triangulation, refine_sequence_value[i]);
      }
      refine_sequence_type.clear();

      make_boundary_face_normals();
      output_boundary_id_areas();

      // Solve forcefield once to get the point inside (not sure if necessary)
      neighborlist->init();

      atom_data->exchange_owned();

      atom_data->exchange_ghost();

      neighborlist->build(true);

      atom_data->reset_owned_acceleration();
      for (auto &&f_custom : force_field_custom)
        f_custom->calculate_acceleration();

      switch (solve_type)
      {
      case 0:
        sg_setup_system();
        sg_assemble_system();
        sg_solve();
        break;

      case 1:
        sa_setup_system();
        sa_assemble_system();
        sa_solve();
        break;

      case 2:
        fg_setup_system();
        fg_assemble_system();
        fg_solve();
        break;

      case 3:
        fa_setup_system();
        fa_assemble_system();
        fa_solve();
        break;
      }
    }

    std::cout << "generate_ml_training_data: Export Boundary Position" << std::endl;
    // Export Boundary Position
    {
      std::ofstream ofs_bound("o_ml_boundary_points");
      // ofs_bound << "X_b"<<","<<"Y_b"<<","<<"Z_b" << "\n";
      for (unsigned int f = 0; f < face_id.size(); ++f)
      { // Potential on the boundary

        if (face_id[f] == 0)
          continue;
        if (std::count(face_id_ignore.begin(), face_id_ignore.end(), face_id[f]) > 0)
          continue;
        ofs_bound << face_center[f][0] << " " << face_center[f][1] << " " << face_center[f][2] << "\n";
      }
    }

    std::cout << "generate_ml_training_data: create grid vector" << std::endl;
    // create all points of the grid vector
    std::vector<Vector<double>> cpoints;
    std::vector<dealii::Point<3>> dpoints;
    auto no_reserve = grid_1d_x->no_points() * grid_1d_y->no_points() * grid_1d_z->no_points();
    cpoints.reserve(no_reserve);
    dpoints.reserve(no_reserve);

    for (unsigned int i = 0; i < grid_1d_x->no_points(); ++i)
    {
      double x = grid_1d_x->give_point(i);
      for (unsigned int j = 0; j < grid_1d_y->no_points(); ++j)
      {
        double y = grid_1d_y->give_point(j);
        for (unsigned int k = 0; k < grid_1d_z->no_points(); ++k)
        {
          double z = grid_1d_z->give_point(k);
          try
          {
            VectorTools::point_value(dof_handler, solution, dealii::Point<3>{x, y, z});
          }
          catch (...)
          {
            // Point is outside of the mesh
            continue;
          }
          cpoints.emplace_back(Vector<double>{x, y, z});
          dpoints.emplace_back(dealii::Point<3>{x, y, z});
        }
      }
    }

    std::cout << "generate_ml_training_data: export ml_inside_points" << std::endl;
    {
      std::ofstream ofs_pt("o_ml_inside_points");
      for (unsigned int i = 0; i < cpoints.size(); ++i)
      {
        ofs_pt << cpoints[i].x << " " << cpoints[i].y << " " << cpoints[i].z << "\n";
      }
    }

    std::cout << "generate_ml_training_data: generate reference data" << std::endl;

    // generate data
    double small_number = 1e-9;
    std::vector<double> boundary_potential_si(face_id.size(), 0);
    std::vector<double> boundary_potential_sm(face_id.size(), 0);
    std::vector<double> boundary_potential_to(face_id.size(), 0);
    auto no_points = cpoints.size();

    std::ofstream ofs_ref_si("o_ml_reference_data_si");
    std::ofstream ofs_ref_sm("o_ml_reference_data_sm");
    std::ofstream ofs_ref_to("o_ml_reference_data_to");

    bool estimate_time = true;
    double t1 = 0, t2 = 0;
    for (unsigned int i = 0; i < no_points; ++i)
    {
      if (estimate_time)
        t1 = get_wall_time();

      // move the charged  particle position
      pos[0] = cpoints[i];

      // calculate force afte the move
      // neighborlist -> init ();

      atom_data->exchange_owned();

      atom_data->exchange_ghost();

      neighborlist->build(true);

      atom_data->reset_owned_acceleration();
      for (auto &&f_custom : force_field_custom)
        f_custom->calculate_acceleration();

      switch (solve_type)
      {
      case 0:
        sg_setup_system();
        sg_assemble_system();
        sg_solve();
        break;

      case 1:
        sa_setup_system();
        sa_assemble_system();
        sa_solve();
        break;

      case 2:
        fg_setup_system();
        fg_assemble_system();
        fg_solve();
        break;

      case 3:
        fa_setup_system();
        fa_assemble_system();
        fa_solve();
        break;
      }

      // Calculate the potential on the boundaries
      for (unsigned int f = 0; f < face_id.size(); ++f)
      { // Potential on the boundary

        if (face_id[f] == 0)
          continue;
        if (std::count(face_id_ignore.begin(), face_id_ignore.end(), face_id[f]) > 0)
          continue;
        // auto p1 = face_center[f];
        double p_sm = 0;
        try
        {
          p_sm = VectorTools::point_value(dof_handler, solution, face_center[f]);
        }
        catch (...)
        {
          error->all(FC_FILE_LINE_FUNC, "face_center point is outside of the mesh.");
          continue;
        }
        double p_si = 0;
        for (auto &&f_custom : force_field_custom)
          p_si += f_custom->potential(caviar::Vector<double>{face_center[f][0], face_center[f][1], face_center[f][2]});

        // ofs_ref << p_si + p_sm << ",";
        // boundary_potential[f] = p_si + p_sm;
        double p_to = p_sm + p_si;
        if (abs(p_sm) < small_number)
          p_sm = 0;
        if (abs(p_si) < small_number)
          p_si = 0;
        if (abs(p_to) < small_number)
          p_to = 0;
        boundary_potential_sm[f] = p_sm;
        boundary_potential_si[f] = p_si;
        boundary_potential_to[f] = p_to;
      }

      // potential at the  points inside
      for (unsigned int j = 0; j < no_points; ++j)
      {

        ofs_ref_sm << cpoints[j].x << " " << cpoints[j].y << " " << cpoints[j].z << " ";
        ofs_ref_si << cpoints[j].x << " " << cpoints[j].y << " " << cpoints[j].z << " ";
        ofs_ref_to << cpoints[j].x << " " << cpoints[j].y << " " << cpoints[j].z << " ";
        // Calculate the potential on the boundaries
        for (unsigned int f = 0; f < face_id.size(); ++f)
        { // Potential on the boundary
          if (face_id[f] == 0)
            continue;
          if (std::count(face_id_ignore.begin(), face_id_ignore.end(), face_id[f]) > 0)
            continue;
          ofs_ref_sm << boundary_potential_sm[f] << " ";
          ofs_ref_si << boundary_potential_si[f] << " ";
          ofs_ref_to << boundary_potential_to[f] << " ";
        }

        double p_sm = 0;
        try
        {
          p_sm = VectorTools::point_value(dof_handler, solution, dpoints[j]);
        }
        catch (...)
        {
          error->all(FC_FILE_LINE_FUNC, "point is outside of the mesh.");
          continue;
        }
        double p_si = 0;
        for (auto &&f_custom : force_field_custom)
          p_si += f_custom->potential(cpoints[j]);

        ofs_ref_sm << p_sm << "\n";
        ofs_ref_si << p_si << "\n";
        ofs_ref_to << p_si + p_sm << "\n";
      }

      // potential at the  points  on the boundaries
      for (unsigned int j = 0; j < face_id.size(); ++j)
      {
        if (face_id[j] == 0)
          continue;
        if (std::count(face_id_ignore.begin(), face_id_ignore.end(), face_id[j]) > 0)
          continue;

        ofs_ref_sm << face_center[j][0] << " " << face_center[j][1] << " " << face_center[j][2] << " ";
        ofs_ref_si << face_center[j][0] << " " << face_center[j][1] << " " << face_center[j][2] << " ";
        ofs_ref_to << face_center[j][0] << " " << face_center[j][1] << " " << face_center[j][2] << " ";

        // Calculate the potential on the boundaries
        for (unsigned int f = 0; f < face_id.size(); ++f)
        { // Potential on the boundary
          if (face_id[f] == 0)
            continue;
          if (std::count(face_id_ignore.begin(), face_id_ignore.end(), face_id[f]) > 0)
            continue;
          ofs_ref_sm << boundary_potential_sm[f] << " ";
          ofs_ref_si << boundary_potential_si[f] << " ";
          ofs_ref_to << boundary_potential_to[f] << " ";
        }

        ofs_ref_sm << boundary_potential_sm[j] << "\n";
        ofs_ref_si << boundary_potential_si[j] << "\n";
        ofs_ref_to << boundary_potential_sm[j] + boundary_potential_si[j] << "\n";
      }

      if (estimate_time)
      {
        t2 = get_wall_time();
        std::cout << "No Points Inside: " << no_points << "\n";
        std::cout << "No rows: " << no_points * no_points << "\n";
        std::cout << "Estimated file size (MB): " << no_points * no_points * (no_points + face_id.size()) * 10 / (1024 * 1024) << "\n";
        std::cout << "Estimate time  (min): " << ((t2 - t1) * no_points) / 60.0 << std::endl;
        estimate_time = false;
      }
    }
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::output_field_vectors(caviar::interpreter::Parser *parser)
  {

    std::string file_name = "o_field_vectors";
    double scale = 1.0;
    double limit = -1.0;
    char field_type = 't';
    bool in_file = true;
    if (in_file == true)
    {
      // removes a warning
    }
    unique::Grid_1D *grid_1d_x = nullptr, *grid_1d_y = nullptr, *grid_1d_z = nullptr;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION

      auto t = token.string_value;

      if (string_cmp(t, "file_name"))
      {
        auto t2 = parser->get_val_token();
        file_name = t2.string_value;
      }
      else if (string_cmp(t, "type"))
      {
        auto t2 = parser->get_val_token();
        if (t2.string_value == "total")
          field_type = 't';
        if (t2.string_value == "smooth")
          field_type = 'm';
        if (t2.string_value == "singular")
          field_type = 'i';
      }
      else if (string_cmp(t, "scale"))
      {
        auto t2 = parser->get_val_token();
        scale = t2.real_value;
      }
      else if (string_cmp(t, "limit"))
      {
        auto t2 = parser->get_val_token();
        limit = t2.real_value;
      }
      else if (string_cmp(t, "grid_1d_x"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_1d_x = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_y"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_1d_y = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_z"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_1d_z = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }

      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    FC_NULLPTR_CHECK(grid_1d_x)
    FC_NULLPTR_CHECK(grid_1d_y)
    FC_NULLPTR_CHECK(grid_1d_z)

    std::ofstream ofs;

    Vector<double> po{0, 0, 0};
    if (position_offset != nullptr)
      po += position_offset->current_value;

    ofs.open(file_name.c_str());

    for (unsigned int i = 0; i < grid_1d_x->no_points(); ++i)
    {
      double x = grid_1d_x->give_point(i);
      for (unsigned int j = 0; j < grid_1d_y->no_points(); ++j)
      {
        double y = grid_1d_y->give_point(j);
        for (unsigned int k = 0; k < grid_1d_z->no_points(); ++k)
        {
          double z = grid_1d_z->give_point(k);

          caviar::Vector<double> field_tot = {0, 0, 0};

          // const dealii::Point<3> r = {x, y, z};
          const dealii::Point<3> r = {x - po.x, y - po.y, z - po.z}; // Is it necessary ??

          dealii::Tensor<1, 3, double> field_sm;

          try
          {
            field_sm = -VectorTools::point_gradient(dof_handler, solution, r);
          }
          catch (...)
          {
            continue;
          }

          if (field_type != 'i')
          {
            field_tot.x += field_sm[0];
            field_tot.y += field_sm[1];
            field_tot.z += field_sm[2];
          }

          if (field_type != 'm')
          {
            const caviar::Vector<double> p{x, y, z};
            caviar::Vector<double> field_si{0, 0, 0};

            for (auto &&f : force_field_custom)
              field_si += f->field(p);

            field_tot += field_si;
          }

          if (limit > 0)
          {
            if (std::abs(field_tot.x) > limit)
              field_tot.x = limit * (field_tot.x / std::abs(field_tot.x));
            if (std::abs(field_tot.y) > limit)
              field_tot.y = limit * (field_tot.y / std::abs(field_tot.y));
            if (std::abs(field_tot.z) > limit)
              field_tot.z = limit * (field_tot.z / std::abs(field_tot.z));
          }

          field_tot = field_tot * scale;

          ofs << x << " " << y << " " << z << " "
              << field_tot.x << " " << field_tot.y << " " << field_tot.z << "\n";
        }
      }
    }

    ofs.close();
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::output_potential_values(caviar::interpreter::Parser *parser)
  {

    std::string file_name = "o_field_vectors";
    // double scale = 1.0;
    // double limit = -1.0;
    char field_type = 't';
    bool in_file = true;
    if (in_file == true)
    {
      // removes a warning
    }
    unique::Grid_1D *grid_1d_x = nullptr, *grid_1d_y = nullptr, *grid_1d_z = nullptr;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION

      auto t = token.string_value;

      if (string_cmp(t, "file_name"))
      {
        auto t2 = parser->get_val_token();
        file_name = t2.string_value;
      }
      else if (string_cmp(t, "type"))
      {
        auto t2 = parser->get_val_token();
        if (t2.string_value == "total")
          field_type = 't';
        if (t2.string_value == "smooth")
          field_type = 'm';
        if (t2.string_value == "singular")
          field_type = 'i';
      }
      else if (string_cmp(t, "grid_1d_x"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_1d_x = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_y"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_1d_y = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_z"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_1d_z = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }

      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    FC_NULLPTR_CHECK(grid_1d_x)
    FC_NULLPTR_CHECK(grid_1d_y)
    FC_NULLPTR_CHECK(grid_1d_z)

    std::string st_ofs_x = file_name + "_x";
    std::string st_ofs_y = file_name + "_y";
    std::string st_ofs_z = file_name + "_z";
    std::string st_ofs_p = file_name + "_p";

    std::ofstream ofs;
    std::ofstream ofs_x, ofs_y, ofs_z, ofs_p;

    Vector<double> po{0, 0, 0};
    if (position_offset != nullptr)
      po += position_offset->current_value;

    ofs.open(file_name.c_str());
    ofs_x.open(st_ofs_x.c_str());
    ofs_y.open(st_ofs_y.c_str());
    ofs_z.open(st_ofs_z.c_str());
    ofs_p.open(st_ofs_p.c_str());

    for (unsigned int i = 0; i < grid_1d_x->no_points(); ++i)
    {
      double x = grid_1d_x->give_point(i);
      for (unsigned int j = 0; j < grid_1d_y->no_points(); ++j)
      {
        double y = grid_1d_y->give_point(j);
        for (unsigned int k = 0; k < grid_1d_z->no_points(); ++k)
        {
          double z = grid_1d_z->give_point(k);

          double potential_tot = 0;

          // const dealii::Point<3> r = {x , y , z };
          const dealii::Point<3> r = {x - po.x, y - po.y, z - po.z}; // Is it necessary ??

          double potential_sm = 0;
          try
          {
            potential_sm = VectorTools::point_value(dof_handler, solution, r);
          }
          catch (...)
          {
            continue;
          }

          if (field_type != 'i')
          {
            potential_tot += potential_sm;
          }

          if (field_type != 'm')
          {
            const caviar::Vector<double> p{x, y, z};
            double potential_si = 0;

            for (auto &&f : force_field_custom)
              potential_si += f->potential(p);

            potential_tot += potential_si;
          }

          ofs << x << " " << y << " " << z << " " << potential_tot << "\n";

          ofs_x << x << "\n";
          ofs_y << y << "\n";
          ofs_z << z << "\n";
          ofs_p << potential_tot << "\n";
        }
      }
    }

    ofs.close();

    ofs_x.close();
    ofs_y.close();
    ofs_z.close();
    ofs_p.close();
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::start_spherical_test()
  {
    FC_NULLPTR_CHECK(test_force_spherical)
    FC_NULLPTR_CHECK(atom_data)
    if (force_field_custom.size() < 1)
      error->all(FC_FILE_LINE_FUNC, "expected an electrostatic forcefield added");
    initialized = false;
    // if (!init_test_force_spherical) {
    //   init_test_force_spherical = true;
    // }
    make_boundary_face_normals();
    output_boundary_id_areas();

    switch (solve_type)
    {
    case 0:
      sg_setup_system();
      sg_assemble_system();
      sg_solve();
      break;

    case 1:
      sa_setup_system();
      sa_assemble_system();
      sa_solve();
      break;

    case 2:
      fg_setup_system();
      fg_assemble_system();
      fg_solve();
      break;

    case 3:
      fa_setup_system();
      fa_assemble_system();
      fa_solve();
      break;
    }
  }
  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::write_spherical_test()
  {

    if (time_step_count % time_step_solve == 0)
    {
      switch (solve_type)
      {
      case 0:
        sg_setup_system();
        sg_assemble_system();
        sg_solve();
        break;

      case 1:
        sa_setup_system();
        sa_assemble_system();
        sa_solve();
        break;

      case 2:
        if (re_do_setup_assemble)
        {
          fg_setup_system();
          fg_assemble_system();
        }
        fg_solve();
        break;

      case 3:
        if (re_do_setup_assemble)
        {
          fa_setup_system();
          fa_assemble_system();
        }
        fa_solve();
        break;
      }
    }

    double c = calculate_induced_charge(0, 1);

    const auto &pos = atom_data->atom_struct_owned.position;
    const auto type_i = atom_data->atom_struct_owned.type[0];
    const auto charge_i = atom_data->atom_type_params.charge[type_i];

    const dealii::Point<3> r = {pos[0].x, pos[0].y, pos[0].z};

    double p_sm = 0.0;
    dealii::Tensor<1, 3, double> f_sm_d;

    // when the point is out of the mesh, we ignore the output
    bool ignore_output = false;

    try
    {
      f_sm_d = -VectorTools::point_gradient(dof_handler, solution, r);
      p_sm = VectorTools::point_value(dof_handler, solution, r);
    }
    catch (...)
    {
      ignore_output = true;
    }

    if (ignore_output)
    {
    }
    else
    {
      //  Vector<double> f_sm = (f_sm_d[0], f_sm_d[1], f_sm_d[2]);
      // auto f_sm_norm = std::sqrt(f_sm*f_sm);
      auto f_sm_norm = std::sqrt(f_sm_d * f_sm_d);

      auto p_an = test_force_spherical->potential(0);
      auto f_an = test_force_spherical->field(0);
      auto f_an_norm = std::sqrt(f_an * f_an);

      auto p_si = force_field_custom[0]->potential(0);
      auto f_si = force_field_custom[0]->field(0);
      auto f_si_norm = std::sqrt(f_si * f_si);

      auto p_err = std::abs((p_si + p_sm - p_an) / p_an);
      auto f_err = std::abs((f_si_norm + f_sm_norm - f_an_norm) / f_an_norm);

      ofs_test_force_spherical << pos[0] << " "
                               << charge_i << " "
                               << c << " "
                               << p_sm << " "
                               << f_sm_norm << " "
                               << p_si << " "
                               << f_si_norm << " "
                               << p_an << " "
                               << f_an_norm << " "
                               << p_si + p_sm << " "
                               << f_si_norm + f_sm_norm << " "
                               << p_err << " "
                               << f_err << "\n";
      ;
    }
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::output_vtk_solution(const std::string filename) const
  {
#if defined(CAVIAR_WITH_MPI)
    if (my_mpi_rank == 0)
    {
#endif

      DataOut<3> data_out;

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "smoothPotential");

      data_out.build_patches();

      std::ofstream output(filename.c_str());
      data_out.write_vtk(output);
#if defined(CAVIAR_WITH_MPI)
    }
#endif
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::output_vtk_solution(const int cycle) const
  {
#if defined(CAVIAR_WITH_MPI)
    if (my_mpi_rank == 0)
    {
#endif
      std::string filename("o_solution-" +
                           Utilities::int_to_string(cycle, zeros_of_mesh_output) +
                           ".vtk");

      DataOut<3> data_out;

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "smoothPotential");

      data_out.build_patches();

      std::ofstream output(filename.c_str());
      data_out.write_vtk(output);
#if defined(CAVIAR_WITH_MPI)
    }
#endif
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::make_boundary_face_normals()
  {
    bool non_planar_face_found = false;
    face_area.clear();
    face_normal.clear();
    face_center.clear();
    face_id.clear();

    std::ofstream ofs;

#if defined(CAVIAR_WITH_MPI)
    if (my_mpi_rank == 0)
#endif
      ofs.open("o_mesh_boundary_normals");

    for (typename Triangulation<3, 3>::active_cell_iterator
             cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
    {
      auto cc = cell->center();

      for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
      {

        if (cell->face(f)->at_boundary())
        {

          auto boundary_id = static_cast<unsigned>(cell->face(f)->boundary_id());

          if (boundary_id_max < boundary_id)
            boundary_id_max = boundary_id;

          auto fc = cell->face(f)->center();

          // n_o: a vector (approximately) to the direction of normal
          auto n_o = cc - fc;

          // this is of type 'dealii::Tensor<1,3> '
          const Tensor<1, 3> v01 = cell->face(f)->vertex(1) - cell->face(f)->vertex(0);
          const Tensor<1, 3> v02 = cell->face(f)->vertex(2) - cell->face(f)->vertex(0);

          Tensor<1, 3> normal = cross_product_3d(v01, v02);

          auto dot = normal * n_o;
          if (dot < 0)
            normal *= -1;
          auto norm_n = std::sqrt(normal * normal);
          normal /= norm_n;

          //

          double area = cell->face(f)->measure();
          if (std::isinf(area) || std::isnan(area))
          {

            const Tensor<1, 3> v03 = cell->face(f)->vertex(3) - cell->face(f)->vertex(0);

            double non_planar_measure = std::abs((v03 * normal) * (v03 * normal) /
                                                 ((v03 * v03) * (v01 * v01) * (v02 * v02)));
            if (non_planar_measure >= 1e-24)
            {
              non_planar_face_found = true;
              // std::cout << "non planar face, measure : "<< non_planar_measure  << "\n";
            }

            const Tensor<1, 3> v12 = cell->face(f)->vertex(2) - cell->face(f)->vertex(1);
            Tensor<1, 3> twice_area = cross_product_3d(v03, v12);
            area = 0.5 * twice_area.norm();
          }
          //

          face_area.push_back(area);
          face_normal.push_back(normal);
          face_center.push_back(fc);
          face_id.push_back(boundary_id);
#if defined(CAVIAR_WITH_MPI)
          if (my_mpi_rank == 0)
#endif
            ofs << fc(0) << " " << fc(1) << " " << fc(2) << " "
                << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
          // std::cout << "fa " << cell->face(f)->measure () << " fi " << boundary_id << "\n";
          // std::cout << "fx:" << cell->face(f)->index() << std::endl;
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

  void Plt_dealii::output_boundary_id_areas()
  {

    std::vector<double> boundary_id_area(boundary_id_max + 1, 0);

    for (unsigned int f = 0; f < face_id.size(); ++f)
    {
      boundary_id_area[face_id[f]] += face_area[f];
      // auto p1 = face_center[f];
    }

    std::ofstream ofs;

#if defined(CAVIAR_WITH_MPI)
    if (my_mpi_rank == 0)
#endif
      ofs.open("o_mesh_boundary_area");

#if defined(CAVIAR_WITH_MPI)
    if (my_mpi_rank == 0)
#endif
      for (unsigned int i = 0; i < boundary_id_area.size(); ++i)
      {
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

  double Plt_dealii::calculate_induced_charge(int t, const int requested_id)
  {

    std::vector<double> induced_charge(boundary_id_max + 1, 0);
    Vector<double> po{0, 0, 0};
    if (position_offset != nullptr)
      po += position_offset->current_value;

    // Not sure if this loop can or should be parallelized with openmp
    // Since there's field() calculation inside, it already has an intrinsic opm
    // parallelized loop inside
    // #ifdef CAVIAR_WITH_OPENMP
    // #pragma omp parallel for
    // #endif
    for (unsigned int f = 0; f < face_id.size(); ++f)
    {

      if (face_id[f] == 0)
        continue;
      if (std::count(face_id_ignore.begin(), face_id_ignore.end(), face_id[f]) > 0)
        continue;
      auto p1 = face_center[f];

      // auto f_sm = - VectorTools::point_gradient (dof_handler, solution, p1);

      dealii::Tensor<1, 3, double> f_sm;
      if (ignore_point_out_of_mesh)
      {
        try
        {
          f_sm = -VectorTools::point_gradient(dof_handler, solution, p1);
        }
        catch (...)
        {
          continue;
        }
      }
      else
      {
        f_sm = -VectorTools::point_gradient(dof_handler, solution, p1);
      }

      caviar::Vector<double> pos_j{p1[0] + po.x, p1[1] + po.y, p1[2] + po.z};

      caviar::Vector<double> f_si{0.0, 0.0, 0.0};
      for (auto &&f_custom : force_field_custom)
        f_si += f_custom->field(pos_j); // OPEN MP IS INSIDE

#ifdef CAVIAR_WITH_MPI
      double f_si_local[3] = {f_si.x, f_si.y, f_si.z};
      double f_si_total[3] = {0, 0, 0};

      MPI_Allreduce(&f_si_local, &f_si_total,
                    3, MPI::DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      auto field_tot = dealii::Point<3>{f_si_total[0], f_si_total[1], f_si_total[2]} + f_sm;
#else
      auto field_tot = dealii::Point<3>{f_si.x, f_si.y, f_si.z} + f_sm;
#endif

      auto field_normal = face_normal[f] * field_tot;

      auto local_q = field_normal * face_area[f]; // this is the charge value times 4*PI*K_el

      // #ifdef CAVIAR_WITH_OPENMP
      // #pragma omp atomic
      // #endif
      induced_charge[face_id[f]] += local_q;
      // std::cout << face_normal[f]<< " " << local_q << " " << face_area[f] << "\n";
    }

    static auto four_pi_k_electrostatic_inv = 1.0 / (4.0 * 3.14159265 * k_electrostatic);

    for (auto &&i : induced_charge)
    {
      i *= four_pi_k_electrostatic_inv; // correcting to the absolute charge value
    }

    ofs_induced_charge << std::setprecision(12); // XXX

    double sum_q = 0; // This is defined here because of CAVIAR_WITH_MPI

#if defined(CAVIAR_WITH_MPI)
    if (my_mpi_rank == 0)
    {
#endif

      if (!induced_charge_id_init)
      {
        induced_charge_id_init = true;
        ofs_induced_charge << "# time ";
        for (unsigned int i = 0; i < induced_charge.size(); ++i)
        {

          // if (induced_charge[i] == 0) continue; // old one
          if (std::count(face_id_ignore.begin(), face_id_ignore.end(), i) > 0)
            continue; // more general

          ofs_induced_charge << " " << i;
        }
        ofs_induced_charge << " "
                           << "sum_q"
                           << " "
                           << "sum_abs_q"
                           << "\n"
                           << std::flush;
      }

      // double sum_q = 0; // This is not defined here because of CAVIAR_WITH_MPI
      double sum_abs_q = 0;
      ofs_induced_charge << t;

      for (unsigned int i = 0; i < induced_charge.size(); ++i)
      {
        if (std::count(face_id_ignore.begin(), face_id_ignore.end(), i) > 0)
          continue;
        auto q = induced_charge[i];
        ofs_induced_charge << " " << q;
        sum_q += q;
        sum_abs_q += std::abs(q);
      }
      ofs_induced_charge << " " << sum_q << " " << sum_abs_q << "\n"
                         << std::flush;

#if defined(CAVIAR_WITH_MPI)
    }
#endif
    if (requested_id != -1)
      return induced_charge[requested_id];
    return sum_q;
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::run()
  {

    if (!initialized)
    {
      initialized = true;

      FC_NULLPTR_CHECK(atom_data)
      domain = atom_data->domain;

      for (unsigned int i = 0; i < refine_sequence_type.size(); ++i)
      {
        if (refine_sequence_type[i] == 0)
          triangulation.refine_global(refine_sequence_value[i]);
        if (refine_sequence_type[i] == 1)
          dealii_functions::refine_boundary(triangulation, refine_sequence_value[i]);
      }
      refine_sequence_type.clear();

      make_boundary_face_normals();
      output_boundary_id_areas();

      switch (solve_type)
      {
      case 0:
        sg_setup_system();
        sg_assemble_system();
        sg_solve();
        break;

      case 1:
        sa_setup_system();
        sa_assemble_system();
        sa_solve();
        break;

      case 2:
        fg_setup_system();
        fg_assemble_system();
        fg_solve();
        break;

      case 3:
        fa_setup_system();
        fa_assemble_system();
        fa_solve();
        break;
      }

      output_vtk_solution(0);

      if (output_induced_charge)
        calculate_induced_charge(0);
    }
    else
    {

      ++time_step_count;

      if (time_step_count % time_step_solve == 0)
      {
        switch (solve_type)
        {
        case 0:
          sg_setup_system();
          sg_assemble_system();
          sg_solve();
          break;

        case 1:
          sa_setup_system();
          sa_assemble_system();
          sa_solve();
          break;

        case 2:
          if (re_do_setup_assemble)
          {
            fg_setup_system();
            fg_assemble_system();
          }
          fg_solve();
          break;

        case 3:
          if (re_do_setup_assemble)
          {
            fa_setup_system();
            fa_assemble_system();
          }
          fa_solve();
          break;
        }
      }

      if (output_vtk && time_step_count % time_step_output_vtk == 0)
        output_vtk_solution(time_step_count);

      if (output_induced_charge && time_step_count % time_step_induced_charge == 0)
        calculate_induced_charge(time_step_count);
    }
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::run_time_profile()
  {

    if (!initialized)
    {
      initialized = true;

      for (unsigned int i = 0; i < refine_sequence_type.size(); ++i)
      {
        if (refine_sequence_type[i] == 0)
          triangulation.refine_global(refine_sequence_value[i]);
        if (refine_sequence_type[i] == 1)
          dealii_functions::refine_boundary(triangulation, refine_sequence_value[i]);
      }
      refine_sequence_type.clear();

      make_boundary_face_normals();
      output_boundary_id_areas();

      double t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0;
      switch (solve_type)
      {
      case 0:
        t1 = get_wall_time();
        sg_setup_system();
        t2 = get_wall_time();
        sg_assemble_system();
        t3 = get_wall_time();
        sg_solve();
        t4 = get_wall_time();
        break;

      case 1:
        t1 = get_wall_time();
        sa_setup_system();
        t2 = get_wall_time();
        sa_assemble_system();
        t3 = get_wall_time();
        sa_solve();
        t4 = get_wall_time();
        break;

      case 2:
        t1 = get_wall_time();
        fg_setup_system();
        t2 = get_wall_time();
        fg_assemble_system();
        t3 = get_wall_time();
        fg_solve();
        t4 = get_wall_time();
        break;

      case 3:
        t1 = get_wall_time();
        fa_setup_system();
        t2 = get_wall_time();
        fa_assemble_system();
        t3 = get_wall_time();
        fa_solve_time_profile();
        t4 = get_wall_time();
        break;
      }

      output_vtk_solution(0);

      t5 = get_wall_time();
      if (output_induced_charge)
        calculate_induced_charge(0);
      t6 = get_wall_time();

      std::cout << "setup: " << (t2 - t1) << "\n";
      std::cout << "assemble: " << (t3 - t2) << "\n";
      std::cout << "solve: " << (t4 - t3) << "\n";
      std::cout << "induced_charge: " << (t6 - t5) << "\n";
    }
    else
    {

      ++time_step_count;
      double t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0;
      if (time_step_count % time_step_solve == 0)
      {

        t1 = get_wall_time();
        t2 = get_wall_time();
        switch (solve_type)
        {
        case 0:
          t1 = get_wall_time();
          sg_setup_system();
          t2 = get_wall_time();
          sg_assemble_system();
          t3 = get_wall_time();
          sg_solve();
          t4 = get_wall_time();
          break;

        case 1:
          t1 = get_wall_time();
          sa_setup_system();
          t2 = get_wall_time();
          sa_assemble_system();
          t3 = get_wall_time();
          sa_solve();
          t4 = get_wall_time();
          break;

        case 2:
          if (re_do_setup_assemble)
          {
            t1 = get_wall_time();
            fg_setup_system();
            t2 = get_wall_time();
            fg_assemble_system();
            t3 = get_wall_time();
          }
          t3 = get_wall_time();
          fg_solve();
          t4 = get_wall_time();
          break;

        case 3:
          if (re_do_setup_assemble)
          {
            t1 = get_wall_time();
            fa_setup_system();
            t2 = get_wall_time();
            fa_assemble_system();
            t3 = get_wall_time();
          }
          t3 = get_wall_time();
          fa_solve_time_profile();
          t4 = get_wall_time();
          break;
        }
      }

      if (output_vtk && time_step_count % time_step_output_vtk == 0)
        output_vtk_solution(time_step_count);

      t5 = get_wall_time();
      if (output_induced_charge && time_step_count % time_step_induced_charge == 0)
        calculate_induced_charge(time_step_count);
      t6 = get_wall_time();

      std::cout << "setup: " << (t2 - t1) << "\n";
      std::cout << "assemble: " << (t3 - t2) << "\n";
      std::cout << "solve: " << (t4 - t3) << "\n";
      std::cout << "induced_charge: " << (t6 - t5) << "\n";
    }
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    domain = atom_data->domain;
    my_mpi_rank = atom_data->get_mpi_rank();
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    if (time_profile)
    {
      double t1 = get_wall_time();

      run_time_profile();
      double t2 = get_wall_time();

      calculate_all_particles_mesh_force_acc();

      double t3 = get_wall_time();

      std::cout << "total_run      : " << (t2 - t1) << "\n";
      std::cout << "mesh_force_acc : " << (t3 - t2) << "\n";
    }
    else
    {
      run();

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
    MPI_Comm_rank (MPI_COMM_WORLD, &me);
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);

    MPI_Barrier (MPI_COMM_WORLD);
    const auto &pos = atom_data -> atom_struct_owned.position;
    const auto &type = atom_data -> atom_struct_owned.type;
    int root = 0;
    while (root < nprocs) {

      unsigned root_pos_size = pos.size();
      MPI_Bcast (&root_pos_size, 1,  MPI::UNSIGNED, root, MPI_COMM_WORLD);
      MPI_Barrier (MPI_COMM_WORLD);

      if (root == me) {


        for (int i = root+1; i < nprocs; ++i) {
          MPI_Send (type.data(), root_pos_size, MPI::UNSIGNED, i, 0, MPI_COMM_WORLD);
          MPI_Send (pos.data(), 3*root_pos_size, MPI::DOUBLE, i, 1, MPI_COMM_WORLD);
        }

        std::vector<Vector<double>> root_acc (root_pos_size,{0,0,0});

        for (int i = root+1; i < nprocs; ++i) {
          MPI_Recv (root_acc.data(), 3*root_pos_size, MPI::DOUBLE,
            root, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (unsigned  i=0;i<pos.size();++i) {
          atom_data -> atom_struct_owned.acceleration [i] += root_acc[i];
        }

      } else {
        if (me > root) {
          std::vector<Vector<double>> root_pos   (root_pos_size);
          std::vector<unsigned> root_type (root_pos_size);

          MPI_Recv (root_type.data(), root_pos_size, MPI::UNSIGNED,
            root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          MPI_Recv (root_pos.data(), 3*root_pos_size, MPI::DOUBLE,
            root, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


          std::vector<Vector<double>> root_acc (root_pos_size,{0,0,0});


          for (unsigned i=0;i<pos.size();++i) {
             const auto type_i = atom_data -> atom_struct_owned.type [i] ;
            const auto mass_inv_i = atom_data -> atom_type_params.mass_inv [ type_i ];
            const auto charge_i = atom_data -> atom_type_params.charge [ type_i ];


            for (unsigned j=i+1;j<root_pos.size();++j) {
               const auto type_j = root_type [j] ;
              const auto mass_inv_j = atom_data -> atom_type_params.mass_inv [ type_j ];
              const auto charge_j = atom_data -> atom_type_params.charge [ type_j ];
              const auto dr = root_pos[j] - pos[i];
              const auto dr_sq = dr*dr;
              const auto dr_norm = std::sqrt(dr_sq);
              const auto force = k_electrostatic * charge_i * charge_j * dr / (dr_sq*dr_norm);
              atom_data -> atom_struct_owned.acceleration [i] -= force * mass_inv_i;
              root_acc [j] += force * mass_inv_j; // no periodic boundary condition yet
            }
          }
        }
      }

      MPI_Barrier (MPI_COMM_WORLD);
      ++root;
    }
    */

#endif
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::calculate_all_particles_mesh_force_acc()
  {
    double virialLocal = 0;
    const auto &pos = atom_data->atom_struct_owned.position;
    caviar::Vector<double> po{0, 0, 0};
    caviar::Vector<int> msd_dummy{0, 0, 0};
    bool bool_dummy = false;

    bool get_pressure_process = atom_data->get_pressure_process();

    if (position_offset != nullptr)
      po += position_offset->current_value;
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < pos.size(); ++i)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      const auto type_i = atom_data->atom_struct_owned.type[i];
      const auto mass_inv_i = atom_data->atom_type_params.mass_inv[type_i];
      const auto charge_i = atom_data->atom_type_params.charge[type_i];

      caviar::Vector<double> p_i = pos[i] - po;
      if (atom_data->atom_struct_owned.molecule_index[i] != -1)
        p_i = domain->fix_position(p_i, msd_dummy, bool_dummy); // bool_dummy and msd_dummy are not used;

      const dealii::Point<3> r = {p_i.x, p_i.y, p_i.z};

      dealii::Tensor<1, 3, double> field;
      if (ignore_point_out_of_mesh)
      {
        try
        {
          field = -VectorTools::point_gradient(dof_handler, solution, r);
        }
        catch (...)
        {
          continue;
        }
      }
      else
      {
        field = -VectorTools::point_gradient(dof_handler, solution, r);
      }
      auto frc = field * charge_i;
      auto force = caviar::Vector<double>{frc[0], frc[1], frc[2]};

      if (get_pressure_process)
      {
        atom_data->add_to_external_virial(force, i);
        // or ???
        // atom_data->add_to_external_virial(force, i, p_i);
      }

      atom_data->atom_struct_owned.acceleration[i] += force * mass_inv_i;

      // if (atom_data->pressure_process)
      //   atom_data->add_to_pressure(force*pos[i]);
    }
    atom_data->virialForce += virialLocal;
  }

  //==================================================
  //==================================================
  //==================================================

  double Plt_dealii::potential(const int)
  {
    error->all(FC_FILE_LINE_FUNC, "Not implemented");
    return 0.0;
  }

  double Plt_dealii::potential(const Vector<double> &v)
  {

    double potential_singular = 0.0;
    for (auto &&f : force_field_custom)
      potential_singular += f->potential(v);

    Vector<double> po{0, 0, 0};
    if (position_offset != nullptr)
      po += position_offset->current_value;

    const dealii::Point<3> r = {v.x - po.x, v.y - po.y, v.z - po.z};

    double potential_smooth = 0.0;
    try
    {
      potential_smooth = VectorTools::point_value(dof_handler, solution, r);
    }
    catch (...)
    {
    }

    return potential_smooth + potential_singular;
  }

  //==================================================
  //==================================================
  //==================================================

  void Plt_dealii::scale_position(double scale_ratio, caviar::Vector<int> scale_axis)
  {
    bool x_axis = (scale_axis.x == 1 ? true : false);
    bool y_axis = (scale_axis.x == 1 ? true : false);
    bool z_axis = (scale_axis.x == 1 ? true : false);

    caviar::Vector<double> scale_ratio_3d{1, 1, 1};
    if (scale_axis.x == 1)
      scale_ratio_3d.x *= scale_ratio;
    if (scale_axis.y == 1)
      scale_ratio_3d.y *= scale_ratio;
    if (scale_axis.z == 1)
      scale_ratio_3d.z *= scale_ratio;

    GridTools::transform(
        [scale_ratio_3d](const Point<3> &in)
        {
          return Point<3>(in[0] * scale_ratio_3d.x, in[1] * scale_ratio_3d.y, in[2] * scale_ratio_3d.z);
        },
        triangulation);

    for (unsigned int i = 0; i < face_center.size(); ++i)
    {
      face_center[i][0] *= scale_ratio_3d.x;

      face_center[i][1] *= scale_ratio_3d.y;

      face_center[i][2] *= scale_ratio_3d.z;
    }
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE

#else

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Plt_dealii::Plt_dealii(CAVIAR *fptr) : caviar::Force_field{fptr}
  {
    error->all(FC_FILE_LINE_FUNC, "Recompile CAVIAR with DEAL.II library.");
  }

  Plt_dealii::~Plt_dealii() {}

  bool Plt_dealii::read(caviar::interpreter::Parser *parser)
  {
    error->all(FC_FILE_LINE_FUNC, "Recompile CAVIAR with DEAL.II library.");
    parser->get_val_token();
    return true;
  }

  void Plt_dealii::calculate_acceleration()
  {
    error->all(FC_FILE_LINE_FUNC, "Recompile CAVIAR with DEAL.II library.");
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
#endif
