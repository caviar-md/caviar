
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

#ifdef CAVIAR_WITH_DEALII

#include "caviar/objects/force_field/plt_dealii.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/unique/time_function.h"

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
#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  /*
  void Plt_dealii::sg_setup_system ()
  {
    dof_handler.distribute_dofs (fe);

  //  std::cout << "   Number of degrees of freedom: "
  //            << dof_handler.n_dofs()
  //            << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
  }



  void Plt_dealii::sg_assemble_system ()
  {
    QGauss<3> quadrature_formula(num_quadrature_points);





    //FEValues<3> fe_values (fe, quadrature_formula,
    //                         update_values | update_gradients |
    //                         update_quadrature_points | update_JxW_values);
    FEValues<3> fe_values (fe, quadrature_formula,
                              update_gradients |  update_JxW_values);




    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    dealii::Vector<double> cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    typename DoFHandler<3>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        cell_matrix = 0;
        cell_rhs = 0;
        for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += (fe_values.shape_grad (i, q_index) *
                                     fe_values.shape_grad (j, q_index) *
                                     fe_values.JxW (q_index));


            }

        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              system_matrix.add (local_dof_indices[i],
                                 local_dof_indices[j],
                                 cell_matrix(i,j));


          }
      }





    std::map<types::global_dof_index,double> boundary_values;


    for (auto&& i : boundary_id_value) {
        VectorTools::interpolate_boundary_values (dof_handler,
          i.first,
          plt_dealii::BoundaryValues(i.second, this),
          boundary_values);
    }

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
  }




  void Plt_dealii::sg_solve ()
  {
  //  SolverControl solver_control (1000, 1e-12*system_rhs.l2_norm());
    SolverControl solver_control (solver_control_maximum_iteration, solver_control_tolerance);
    SolverCG<> solver (solver_control);
  //  SolverBicgstab<> solver (solver_control);
    solver.solve (system_matrix, solution, system_rhs,
                  PreconditionIdentity());



  //  std::cout << "   " << solver_control.last_step()
  //            << " CG iterations needed to obtain convergence."
  //            << std::endl;
  }
  */

  void Plt_dealii::sg_setup_system()
  {
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }

  void Plt_dealii::sg_assemble_system()
  {
    QGauss<3> quadrature_formula(num_quadrature_points);

    // const RightHandSide<3> right_hand_side;

    // FEValues<3> fe_values (fe, quadrature_formula,
    //                          update_values   | update_gradients |
    //                          update_quadrature_points | update_JxW_values);

    FEValues<3> fe_values(fe, quadrature_formula,
                          update_gradients | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    dealii::Vector<double> cell_rhs(dofs_per_cell);
    cell_rhs = 0;

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<3>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();

    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      cell_matrix = 0;

      for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            cell_matrix(i, j) += (fe_values.shape_grad(i, q_index) *
                                  fe_values.shape_grad(j, q_index) *
                                  fe_values.JxW(q_index));
        }

      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));
      }
    }

    std::map<types::global_dof_index, double> boundary_values;

    for (auto &&i : boundary_id_value)
    {
      VectorTools::interpolate_boundary_values(dof_handler,
                                               i.first,
                                               plt_dealii::BoundaryValues(i.second, this),
                                               boundary_values);
    }

    for (auto &&i : boundary_id_time_function)
    {
      auto fvalue = i.second->value();
      VectorTools::interpolate_boundary_values(dof_handler,
                                               i.first,
                                               plt_dealii::BoundaryValues(fvalue, this),
                                               boundary_values);
    }

    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
  }

  void Plt_dealii::sg_solve()
  {
    SolverControl solver_control(solver_control_maximum_iteration, solver_control_tolerance);
    SolverCG<> solver(solver_control);

    if (use_preconditioner)
    {
      PreconditionSSOR<> preconditioner;
      preconditioner.initialize(system_matrix, 1.2);

      solver.solve(system_matrix, solution, system_rhs,
                   preconditioner);
    }
    else
    {
      solver.solve(system_matrix, solution, system_rhs,
                   PreconditionIdentity());
    }
  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
#endif
