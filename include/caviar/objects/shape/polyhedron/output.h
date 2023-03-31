
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

#ifndef CAVIAR_OBJECTS_SHAPE_POLYHEDRON_OUTPUT_H
#define CAVIAR_OBJECTS_SHAPE_POLYHEDRON_OUTPUT_H

#include "caviar/utility/pointers.h"

CAVIAR_NAMESPACE_OPEN

namespace shape {
namespace polyhedron {
struct Polyhedron;
class Output : public Pointers{
public:
  Output (class CAVIAR *);
  ~Output ();
  
  void mesh_povray (const shape::polyhedron::Polyhedron &, std::string file = "o_mesh.pov"); // povray output mesh ".pov"
  void mesh_tcl (const shape::polyhedron::Polyhedron &, std::string file = "o_mesh.tcl"); // vfptr output mesh ".tcl"
  void normals_tcl (const shape::polyhedron::Polyhedron &, std::string file = "o_normals.tcl" ); // vfptr output normals ".tcl"
  void normals_vectors (const shape::polyhedron::Polyhedron &, std::string file = "o_normals.txt" ); // vfptr output normals ".txt"
  void edges_tcl (const shape::polyhedron::Polyhedron &, std::string file = "o_edges.tcl"); // vfptr output edges  ".tcl"
   
  //void mesh_povray (const std::vector<polyhedron::Polyhedron> &); // povray output mesh ".pov"
//  void mesh_vfptr (const std::vector<polyhedron::Polyhedron> &); // vfptr output mesh ".tcl"
//  void normals_vfptr (const std::vector<polyhedron::Polyhedron> &); // vfptr output normals ".tcl"
  //void edges_vfptr (const std::vector<polyhedron::Polyhedron> &); // vfptr output edges  ".tcl"
   

};
} //polyhedron
} //shape


CAVIAR_NAMESPACE_CLOSE

#endif
