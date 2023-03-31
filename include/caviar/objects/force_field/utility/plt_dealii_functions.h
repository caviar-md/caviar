
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_PLTDEALIIFUNCTIONS_H
#define CAVIAR_OBJECTS_FORCEFIELD_PLTDEALIIFUNCTIONS_H


// XXX: Note that (for now!) this file should be included after Deal.II usual 
// inclusions.

namespace caviar {

namespace force_field {
namespace dealii_functions {

//===================================================
//===================================================
//===================================================

template <typename T>
void rotate_and_add (T &triangulation, const double angle, const int axis, const double merge_toll) {
  T tria2;
  tria2.copy_triangulation (triangulation);
  GridTools::rotate (angle, axis, tria2);
  match_boundary_vertices (tria2, merge_toll);     
  GridGenerator::merge_triangulations (triangulation, tria2, triangulation);        
}

//===================================================
//===================================================
//===================================================

template <typename T>
void rotate_and_add_reserve(T &triangulation, T &tria_reserve, const double angle, const int axis, const double merge_toll) {
    T tria2;
    tria2.copy_triangulation (tria_reserve);
    GridTools::rotate (angle, axis, tria2);
    match_boundary_vertices (tria2, merge_toll);     
    GridGenerator::merge_triangulations (triangulation, tria2, triangulation);        
}

//===================================================
//===================================================
//===================================================

template <typename T>
void match_boundary_vertices (T &triangulation, T &tria2,
                                               const double merge_toll) {
  int no_matched_vertices = 0;
  int no_corrected_vertices = 0;
    
  for (typename Triangulation<3,3>::active_cell_iterator
      cell2=tria2.begin_active(); cell2!=tria2.end(); ++cell2) {

    for (unsigned int f2=0; f2<GeometryInfo<3>::faces_per_cell; ++f2) {

      if (cell2->face(f2)->at_boundary()) {
          
        for (typename Triangulation<3,3>::active_cell_iterator
            cell=triangulation.begin_active();cell!=triangulation.end();++cell) { 

          for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f) {
            if (cell->face(f)->at_boundary()) {
      
              for (unsigned int i=0; i<4; ++i)  { 
                const Point<3> p1 = cell->face(f)->vertex(i);
                bool point_match_found = false;                             
                for (unsigned int j=0; j<4; ++j)  {           
                  const Point<3> p2 = cell2->face(f2)->vertex(j);
                  const double distance = p2.distance(p1);
                    
                  if (distance < merge_toll) {
                    ++no_matched_vertices;
                    point_match_found = true;
                    if (distance > 0) {
                      ++no_corrected_vertices;
                      cell2->face(f2)->vertex(j) =  cell->face(f)->vertex(i);
                    }
                    break;                        
                  }
                }

                if (!point_match_found) break;
              }
            }
          }
        }           
      }        
    }
  }
    
}




//===================================================
//===================================================
//===================================================

template <typename T>
bool dealii_grid_generator_hyper_rectangle(class CAVIAR *fptr, interpreter::Parser *parser, T &triangulation) {
  auto object_container = fptr->object_container;
  auto error = fptr->error;
  double xmin = 0, xmax = 1 , ymin = 0, ymax = 1, zmin = 0, zmax = 1;
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"xmin") ) {
      GET_OR_CHOOSE_A_REAL(xmin,"","")
    } else if (string_cmp(t,"xmax") ) {
      GET_OR_CHOOSE_A_REAL(xmax,"","")
    } else if (string_cmp(t,"ymin") ) {
      GET_OR_CHOOSE_A_REAL(ymin,"","")
    } else if (string_cmp(t,"ymax") ) {
      GET_OR_CHOOSE_A_REAL(ymax,"","")
    } else if (string_cmp(t,"zmin") ) {
      GET_OR_CHOOSE_A_REAL(zmin,"","")
    } else if (string_cmp(t,"zmax") ) {
      GET_OR_CHOOSE_A_REAL(zmax,"","")
    } 
    else FC_ERR_UNDEFINED_VAR(t)
  }

  if (xmin > xmax) error->all (FC_FILE_LINE_FUNC, "xmin > xmax");
  if (ymin > ymax) error->all (FC_FILE_LINE_FUNC, "ymin > ymax");
  if (zmin > zmax) error->all (FC_FILE_LINE_FUNC, "zmin > zmax");

  const dealii::Point< 3 >  p1 (xmin, ymin, zmin);
  const dealii::Point< 3 >  p2 (xmax, ymax, zmax);

  GridGenerator::hyper_rectangle 	(triangulation, p1, p2, true);

  return in_file;
}


  
//===================================================
//===================================================
//===================================================

template <typename T>
void set_spherical_manifold(T &triangulation) {
  const Point<3> center (0,0,0);
  static const SphericalManifold<3> manifold_description(center);
  triangulation.set_manifold (1, manifold_description);
  triangulation.set_all_manifold_ids(1);
}

//===================================================
//===================================================
//===================================================
 
template <typename T>
bool dealii_grid_generator_hyper_ball (class CAVIAR *fptr, interpreter::Parser *parser, T &triangulation)
{

  auto object_container = fptr->object_container;
  auto error = fptr->error;
  double radius = 1.0;
  caviar::Vector<double> center {0,0,0};
  bool in_file = true;

  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"center") ) {
      GET_OR_CHOOSE_A_REAL_3D_VECTOR(center, "", "");
    } else if (string_cmp(t,"radius") ) {
      GET_OR_CHOOSE_A_REAL(radius,"","")
    } else FC_ERR_UNDEFINED_VAR(t)
  }

  if (radius < 0) error->all (FC_FILE_LINE_FUNC, "radius < 0");



  const Point<3> p_center (center.x, center.y, center.z); 

  GridGenerator::hyper_ball   (triangulation, p_center,    radius);
 
  static const SphericalManifold<3> manifold_description(p_center);
  triangulation.set_manifold (2, manifold_description);
  //triangulation.set_all_manifold_ids(2);

    for (typename Triangulation<3,3>::active_cell_iterator
         cell=triangulation.begin_active(); cell!=triangulation.end(); ++cell) {
      
      for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f) {
        
        if (cell->face(f)->at_boundary()) {
          cell->set_all_manifold_ids(2);
        }
      }
    }
/*
  for (typename Triangulation<3,3>::active_cell_iterator
         cell=triangulation.begin_active(); cell!=triangulation.end(); ++cell) {
    cell->set_all_boundary_ids (0);
  }
*/
/*
  std::cout << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: "
            << triangulation.n_cells()
            << std::endl;
*/

    for (typename Triangulation<3,3>::active_cell_iterator
         cell=triangulation.begin_active(); cell!=triangulation.end(); ++cell) {
      
      for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f) {
        
        if (cell->face(f)->at_boundary()) {
          cell->face(f)->set_boundary_id (1);
          //cell->set_all_boundary_ids (1);
        }
      }
    }

  return in_file;
}


//===================================================
//===================================================
//===================================================

template <typename T>
bool dealii_grid_generator(class CAVIAR *fptr, interpreter::Parser *parser, T &triangulation) {
  auto error = fptr->error;
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"hyper_rectangle")) {
      dealii_grid_generator_hyper_rectangle(fptr, parser, triangulation);
      return in_file;
    } else if (string_cmp(t,"hyper_ball")) {
      dealii_grid_generator_hyper_ball(fptr, parser, triangulation);
      return in_file;
    } else FC_ERR_UNDEFINED_VAR(t)
  }
  return in_file;
}

//==================================================
//==================================================
//==================================================

/*
static
void set_boundary ()
{
// To be added at UNV_HANDLER class if possible
// 

  const double sn = 0.001;
      
  Triangulation<3>::active_cell_iterator cell = triangulation.begin_active();
  Triangulation<3>::active_cell_iterator endc = triangulation.end();    
  for (; cell!=endc; ++cell) {    
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
      if (cell->face(f)->at_boundary()) {

      }          
    }
  }

}
*/

//==================================================
//==================================================
//==================================================

template <typename T>
void refine_boundary (T &triangulation, const unsigned int refine_levels)
{

  for (unsigned int j = 1; j <= refine_levels; ++j) {
    std::cout << "refine level " << j << std::endl;
    Triangulation<3>::active_cell_iterator cell = triangulation.begin_active();
    Triangulation<3>::active_cell_iterator endc = triangulation.end();
    
    for (; cell!=endc; ++cell) {    
      for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary()) {
          if (cell->face(f)->boundary_id()!=0)
            cell->set_refine_flag ();          
        }          
    }
    triangulation.execute_coarsening_and_refinement ();   
  }
}

//===================================================
//===================================================
//===================================================

template <typename T>
void read_domain(T &triangulation,T &tria_reserve, std::vector<std::string> &unv_mesh_filename, const double merge_toll=1e-4)  {

  if (unv_mesh_filename.size()==0) {
    //error->all (FC_FILE_LINE_FUNC_PARSE, "unv_mesh_filename is empty. Add a file to the list.");
    std::cout << "unv_mesh_filename is empty. Add a file to the list." << std::endl;
    return;
  }
  int tot_no_matched=0; 
  int tot_no_corrected=0;


  std::ifstream in;

  in.open(unv_mesh_filename[0].c_str());

  GridIn<3,3> gi;
  gi.attach_triangulation(triangulation);
  gi.read_unv(in);


  std::cout << "read_domain: " << unv_mesh_filename[0] << std::endl;
  unsigned int tot_mesh_num = unv_mesh_filename.size();
  for (unsigned int i = 1; i <  tot_mesh_num; ++i) { // BUG
    std::cout << "read_domain: " << unv_mesh_filename[i] << std::endl;

      tria_reserve.clear();
     
      std::ifstream in;
      
      in.open(unv_mesh_filename[i].c_str());

      GridIn<3,3> gi;
      gi.attach_triangulation(tria_reserve); 
      gi.read_unv(in);
      
      match_boundary_vertices (triangulation, tria_reserve, merge_toll);  
       
      GridGenerator::merge_triangulations (triangulation, tria_reserve, triangulation);
    
  }
    

  std::cout << "merge_tollerance_of_vertices:       " << merge_toll << std::endl;
  std::cout << "total_no_matched_vertices:   " << tot_no_matched << std::endl;
  std::cout << "total_no_corrected_vertices: " << tot_no_corrected << std::endl;

  unv_mesh_filename.clear();

} 

//==================================================
//==================================================
//==================================================

template <typename T_tria, typename T_sol, typename T_dof>
void refine_grid_adaptive (T_tria &triangulation,
    T_sol &solution, T_dof &dof_handler, const double param1=0.3,
    const double param2=0.03) {

  dealii::Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  KellyErrorEstimator<3>::estimate (dof_handler,
                                      QGauss<3-1>(3),
                                      typename FunctionMap<3>::type(),
                                      solution,
                                      estimated_error_per_cell);
  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   param1, param2);
  triangulation.execute_coarsening_and_refinement ();  
  
}

//==================================================
//==================================================
//==================================================

} //dealii_functions
} //force_field

} // namespace caviar
#endif
