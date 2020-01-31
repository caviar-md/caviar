
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

#ifndef UDN_2467_H
#define UDN_2467_H

#include <vector>

namespace mesh_modifier {

struct Universal_dataset_number_2467 {
/*
Record 1:        FORMAT(8I10)
                 Field 1       -- group number
                 Field 2       -- active constraint set no. for group
                 Field 3       -- active restraint set no. for group
                 Field 4       -- active load set no. for group
                 Field 5       -- active dof set no. for group
                 Field 6       -- active temperature set no. for group
                 Field 7       -- active contact set no. for group
                 Field 8       -- number of entities in group

Record 2:        FORMAT(20A2)
                 Field 1       -- group name

Record 3-N:      FORMAT(8I10)
                 Field 1       -- entity type code
                 Field 2       -- entity tag
                 Field 3       -- entity node leaf id.
                 Field 4       -- entity component/ ham id.
                 Field 5       -- entity type code
                 Field 6       -- entity tag
                 Field 7       -- entity node leaf id.
                 Field 8       -- entity component/ ham id.

Repeat record 3 for all entities as defined by record 1, field 8.
Records 1 thru n are repeated for each group in the model.
Entity node leaf id. and the component/ ham id. are zero for all
entities except "reference point", "reference point series"
and "coordinate system".
*/
 std::vector<unsigned int> record1;
 unsigned int record2;
 std::vector<unsigned int> record3;
};


}

#endif
/*

          Permanent group entity type codes

    Entity Type Code        Entity Description

           1                coordinate system
           2                data surface thickness
           3                force on point
           4                force on edge
           5                traction on face
           6                pressure on face
           7                nodes
           8                finite elements
           9                dof sets, dof entities
          10                constraint sets, coupled dofs
          11                constraint sets, mpc equations
          12                restraint sets, nodal displacements
          13                restraint sets, nodal temperatures
          14                load sets, nodal forces
          15                load sets, nodal temperatures
          16                load sets, nodal heat sources/sinks
          17                load sets, face pressures
          18                load sets, edge pressures
          19                load sets, face heat fluxes
          20                load sets, edge heat fluxes
          21                load sets, face heat convections
          22                load sets, edge heat convections
          23                load sets, face heat radiations
          24                load sets, edge heat radiations
          25                load sets, element heat generations
          26                load sets, beam temperatures
          27                trace lines
          28                beam force
          29                beam distributed load
          30                data surface
          31                data curve
          32                displacement on point (restraint)
          33                displacement on edge (restraint)
          34                displacement on surface (restraint)
          35                temperature on point (restraint)
          36                temperature on edge (restraint) 
          37                temperature on face (restraint)
          38                temperature on point (temperature)
          39                temperature on edge (temperature)
          40                temperature on face (temperature)
          41                heat source on point
          42                heat flux on edge
          43                convection on edge
          44                radiation on edge
          45                heat flux on face
          46                convection on face
          47                radiation on face
          48                geometry contact region
          49                fe contact region
          50                contact pair
          51                kinematic dof on point
          52                kinematic dof on edge
          53                kinematic dof on face
          54                element definition
          55                anchor node
          56                edge dependancy mesh definition
          57                fem point connector
          58                fem area connector
          59                vertex
          60                edge
          61                face
          62                region
          63                wireframe connector
          64                wireframe curve
          65                wireframe section
          66                wireframe region
          67                reference point
          68                reference point series
          69                centerpoint

Example:

    -1
  2452
         1         0         1         1         0         5
PERMANENT GROUP1
         7        25         0         0         8         1         0         0
        12         1         0         0        12         2         0         0
        14        22         0         0
    -1
*/

