
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

#ifndef CAVIAR_OBJECTS_SHAPE_POLYHEDRON_POLYHEDRON_H
#define CAVIAR_OBJECTS_SHAPE_POLYHEDRON_POLYHEDRON_H

#include "caviar/utility/caviar_config.h"
#include "caviar/utility/vector.h"
#include "caviar/utility/types.h"
#include <vector>
#include <map>

namespace caviar {

namespace shape {
namespace polyhedron {
struct Polyhedron {

  Polyhedron () : nx_part{1}, ny_part{1}, nz_part{1}, grid_tol{0.0}, thickness{1.0} {}

  /**
   * contains cartesian coordinates
   */
  std::vector<Vector<Real_t>> vertex;  

  /**
   * made in "merge_vertices()" it's a map to the index of possible similar
   * vertex with the lower index
   */
  std::vector<std::vector<unsigned int>> vertex_map;

  /**
   * contains vertex indices of ngons
   */
  std::vector<std::vector<unsigned int>> face;


  /**
   * contains indices of the faces. It can be imported from UNV files or set...
   * The number '-1' is the default and invalid face_id.
   */
  std::vector<int> face_id;


  /**
   * normal of faces. The direction has to be set by the user
   */
  std::vector<Vector<Real_t>> normal;

  /**
   *
   * addition of normals of neighborlist faces for each edge
   *  Locally Defined 
   * std::vector<std::vector<Vector<Real_t>>> edge_norms1;
   *
   * face[][j]-face[][j+1] 
   *  Locally Defined 
   * std::vector<std::vector<Vector<Real_t>>> edge_norms2;
   * 
   * edge_norms1 cross edge_norms2
   */
  std::vector<std::vector<Vector<Real_t>>> edge_norms3; 

  /**
   * first: the edge that has vertex[i],vertex[j] ... second: face[m],face[n] 
   * that has this edge
   */
  std::map<std::vector<unsigned int>,std::vector<unsigned int>> edges; 

  /**
   * contains face indices within the grid;
   */
  std::vector<std::vector<std::vector<std::vector<unsigned int>>>> grid; 

  /**
   * highest and lowest coordinates of vertices; used in grid;
   */
  Real_t xlo, ylo, zlo, xhi, yhi, zhi; 

  /**
   * used in grid;
   */
  Real_t dx_part, dy_part, dz_part; 

  /**
   * used in grid;
   */
  unsigned int nx_part, ny_part, nz_part;   

  /**
   * a tolerance value used in grid;
   */
  Real_t grid_tol;

  /*  
   * it defines the maximum thickness of each surface. If a particle go farther
   * than it inside, it won't be inside anymore.        
   */
  Real_t thickness; 


};
} //polyhedron
} //shape


} // namespace caviar

#endif
