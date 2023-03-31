
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

#include "mesh_modifier.h"

int main()
{

  // example e1
  mesh_modifier::Mesh_modifier mesh_modifier;

  mesh_modifier.import("../examples/e1/mesh-raw/Mesh_1.unv", true);

  mesh_modifier.remove_hexa_internals(1e-8);

  std::vector<Point_condition> pc_e1, pc_e2;

  pc_e1.push_back(Point_condition("x", "<", 0.01));
  pc_e1.push_back(Point_condition("z", ">", 0.01));
  pc_e1.push_back(Point_condition("z", "<", 4.99));

  pc_e2.push_back(Point_condition("x", ">", 19.99));

  mesh_modifier.add_face_to_group_with_condition("2", &pc_e1);
  mesh_modifier.add_face_to_group_with_condition("5", &pc_e2);

  mesh_modifier.export_("../examples/e1/mesh-fixed/Mesh_1_fixed.unv", false);

  /*
    // here we scale some points
    std::vector<Point_condition> pc_s;

    // vertex selection
    pc_s.push_back( Point_condition("x",">",28.9) );
    pc_s.push_back( Point_condition("x","<",29.1) );

    // scale center
    Vector<double> v_center {29,0,0};

    // scale vector
    Vector<double> v_scale {1.0, 1.5/2.0,  1.5/2.0};
    mesh_modifier.scale_vertices(&pc_s, v_center, v_scale);

  */

  /*
    // here we export the boundaries created from UNV mesh
    //mesh_modifier.export_vtk_boundary("../mesh/mesh_raw_pore_2seg_fixed.vtk");
    //mesh_modifier.export_stl_boundary("../mesh/mesh_raw_pore_2seg_fixed.stl");
  */

  /*
    // Here we import a tetra mesh with internal faces and change it to a
    // hexahedral mesh
    mesh_modifier.import("../mesh/Mesh_box_gr.unv", true);
    mesh_modifier.convert_tetra_to_hexa ();
    mesh_modifier.remove_hexa_internals (1e-9);
    mesh_modifier.export_("../mesh/test.unv", false);
  */

  /*
    // adding meshes and removing internal faces, then making a VTK file
    mesh_modifier.import ("../mesh/Mesh_5_merge_g.unv", true);
    mesh_modifier.import ("../mesh/Mesh_6_merge_g.unv", true);
    mesh_modifier.import ("../mesh/Mesh_7_merge_g.unv", true);
    mesh_modifier.import ("../mesh/Mesh_8_merge_g.unv", true);
    mesh_modifier.remove_hexa_internals (1e-9);
    mesh_modifier.export_vtk("../mesh/test.vtk");
  */

  /*
    // import multiple tetra, merge them, convert to hexa, remove internal faces
    mesh_modifier.import("../mesh/Mesh_4_compound_tetra.unv", true);
    mesh_modifier.import("../mesh/Mesh_1.unv", true);
    mesh_modifier.import("../mesh/Mesh_2.unv", true);
    mesh_modifier.merge_all (true);
    mesh_modifier.convert_tetra_to_hexa ();
    mesh_modifier.remove_hexa_internals (1e-9);
    mesh_modifier.export_("../mesh/test.unv", false);
  */
}
