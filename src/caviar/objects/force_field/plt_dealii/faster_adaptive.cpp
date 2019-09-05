
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

#include <cmath>
#include <iomanip>

namespace caviar {
namespace objects {
namespace force_field {



void Plt_dealii::fa_setup_system ()
{

  dof_handler.distribute_dofs (fe);
  
  /*
  std::cout << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;
  */
                    
  solution.reinit (dof_handler.n_dofs()); 
  
  system_rhs.reinit (dof_handler.n_dofs());

  constraints.clear ();
  
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);

  constraints.close ();  

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs = */ false);
  
  sparsity_pattern.copy_from(dsp);
  
  system_matrix.reinit (sparsity_pattern);
}



void Plt_dealii::fa_assemble_system ()
{
  const QGauss<3> quadrature_formula(num_quadrature_points);


  //FEValues<3> fe_values (fe, quadrature_formula,
  //                         update_values | update_gradients |
  //                         update_quadrature_points | update_JxW_values);


  FEValues<3> fe_values (fe, quadrature_formula,
                            update_gradients |  update_JxW_values);
                           
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
  dealii::Vector<double> cell_rhs (dofs_per_cell);
  cell_rhs = 0;
      
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  typename DoFHandler<3>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
 
      fe_values.reinit (cell);
      cell_matrix = 0;

      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad (i, q_index) *
                                   fe_values.shape_grad (j, q_index) *
                                   fe_values.JxW (q_index));


          }


      cell->get_dof_indices (local_dof_indices);

      constraints.distribute_local_to_global (cell_matrix,
                                              cell_rhs,
                                              local_dof_indices,
                                              system_matrix,
                                              system_rhs);
                                              
    }

                                                  
}

//==================================================
//==================================================
//==================================================


void Plt_dealii::fa_solve () // XXX note of two of them
{

  std::map<types::global_dof_index,double> boundary_values;
  
  for (auto&& i : boundary_id_value) {
      VectorTools::interpolate_boundary_values (dof_handler,
        i.first,
        plt_dealii::BoundaryValues(i.second, this),
        boundary_values); 
  }                                            
                                            

  // matrix and boundary value constraints
  FilteredMatrix<dealii::Vector<double> > filtered_A (system_matrix);
  filtered_A.add_constraints (boundary_values);


  // set up a linear solver
  SolverControl solver_control (solver_control_maximum_iteration, solver_control_tolerance, false, false);


  //SolverControl solver_control (1000, 1.e-10);  
  GrowingVectorMemory<dealii::Vector<double> > mem;
  SolverCG<dealii::Vector<double> > solver (solver_control, mem);

 
  
  if (use_preconditioner) {  
    // set up a preconditioner object
    PreconditionJacobi<SparseMatrix<double> > prec;
    prec.initialize (system_matrix, preconditioner_relaxation);
    FilteredMatrix<dealii::Vector<double> > filtered_prec (prec);
    filtered_prec.add_constraints (boundary_values);
  
  
    // compute modification of right hand side
    auto system_rhs_tmp = system_rhs;
    filtered_A.apply_constraints (system_rhs_tmp);

    // solve for solution vector x
    solver.solve (filtered_A, solution, system_rhs_tmp, filtered_prec);
  } else {
    error->all(FC_FILE_LINE_FUNC,"not implemented yet. Use it with preconditioner.");
  }


  constraints.distribute (solution);

}


//==================================================
//==================================================
//==================================================


void Plt_dealii::fa_solve_time_profile () // XXX note the time profile
{

    clock_t t1=0, t2=0, t3=0, t4=0, t5=0, t6=0;
    clock_t t7=0, t8=0, t9=0, t10=0, t11=0, t12=0;
    clock_t t13=0, t14=0, t15=0, t16=0, t17=0;
  std::map<types::global_dof_index,double> boundary_values;
  t1 = clock();  
  for (auto&& i : boundary_id_value) {
          t2 = clock();
      VectorTools::interpolate_boundary_values (dof_handler,
        i.first,
        plt_dealii::BoundaryValues(i.second, this),
        boundary_values); 
          t3 = clock();
  std::cout << "boundary_id " << i.second << "  : "    << (double)(t3 - t2)/CLOCKS_PER_SEC << "\n";
  }                                            
  t4 = clock();                                            
  std::cout << "total_boundary_id: "   << (double)(t4 - t1)/CLOCKS_PER_SEC << "\n";

  // matrix and boundary value constraints
  FilteredMatrix<dealii::Vector<double> > filtered_A (system_matrix);
          t5 = clock();
  filtered_A.add_constraints (boundary_values);
          t6 = clock();

  // set up a linear solver
  SolverControl solver_control (solver_control_maximum_iteration, solver_control_tolerance, false, false);
          t7 = clock();

  //SolverControl solver_control (1000, 1.e-10);  
  GrowingVectorMemory<dealii::Vector<double> > mem;
          t8 = clock();
  SolverCG<dealii::Vector<double> > solver (solver_control, mem);
          t9 = clock();

 
  
  if (use_preconditioner) {  
    // set up a preconditioner object
          t10 = clock();
    PreconditionJacobi<SparseMatrix<double> > prec;
    prec.initialize (system_matrix, preconditioner_relaxation);
          t11 = clock();
    FilteredMatrix<dealii::Vector<double> > filtered_prec (prec);
          t12 = clock();
    filtered_prec.add_constraints (boundary_values);
            t13 = clock();
  
    // compute modification of right hand side
    auto system_rhs_tmp = system_rhs;
          t14 = clock();
    filtered_A.apply_constraints (system_rhs_tmp);
          t15 = clock();
    // solve for solution vector x
    solver.solve (filtered_A, solution, system_rhs_tmp, filtered_prec);
          t16 = clock();
  } else {
    error->all(FC_FILE_LINE_FUNC,"not implemented yet. Use it with preconditioner.");
  }

 
  constraints.distribute (solution);
          t17 = clock();
  std::cout << "t5  - t4: "   << (double)(t5 - t4)/CLOCKS_PER_SEC << "\n";
  std::cout << "t6  - t5: "   << (double)(t6 - t5)/CLOCKS_PER_SEC << "\n";
  std::cout << "t7  - t6: "   << (double)(t7 - t6)/CLOCKS_PER_SEC << "\n";
  std::cout << "t8  - t7: "   << (double)(t8 - t7)/CLOCKS_PER_SEC << "\n";
  std::cout << "t9  - t8: "   << (double)(t9 - t8)/CLOCKS_PER_SEC << "\n";
  std::cout << "t10 - t9: "   << (double)(t10 - t9)/CLOCKS_PER_SEC << "\n";
  std::cout << "t11 - t10: "   << (double)(t11 - t10)/CLOCKS_PER_SEC << "\n";
  std::cout << "t12 - t11: "   << (double)(t12 - t11)/CLOCKS_PER_SEC << "\n";
  std::cout << "t13 - t12: "   << (double)(t13 - t12)/CLOCKS_PER_SEC << "\n";
  std::cout << "t14 - t13: "   << (double)(t14 - t13)/CLOCKS_PER_SEC << "\n";
  std::cout << "t15 - t14: "   << (double)(t15 - t14)/CLOCKS_PER_SEC << "\n";
  std::cout << "t16 - t15: "   << (double)(t16 - t15)/CLOCKS_PER_SEC << "\n";
  std::cout << "t17 - t16: "   << (double)(t17 - t16)/CLOCKS_PER_SEC << "\n";

}


} //force_field
} //objects
} // namespace caviar
#endif
