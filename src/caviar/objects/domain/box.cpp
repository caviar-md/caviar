
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

#include "caviar/objects/domain/box.h"
#include "caviar/interpreter/communicator.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace domain
{

  Box::Box(CAVIAR *fptr) : Domain{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
  }

  bool Box::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "lower_global.x") || string_cmp(t, "xmin"))
      {
        GET_OR_CHOOSE_A_REAL(lower_global.x, "", "")
      }
      else if (string_cmp(t, "upper_global.x") || string_cmp(t, "xmax"))
      {
        GET_OR_CHOOSE_A_REAL(upper_global.x, "", "")
      }
      else if (string_cmp(t, "lower_global.y") || string_cmp(t, "ymin"))
      {
        GET_OR_CHOOSE_A_REAL(lower_global.y, "", "")
      }
      else if (string_cmp(t, "upper_global.y") || string_cmp(t, "ymax"))
      {
        GET_OR_CHOOSE_A_REAL(upper_global.y, "", "")
      }
      else if (string_cmp(t, "lower_global.z") || string_cmp(t, "zmin"))
      {
        GET_OR_CHOOSE_A_REAL(lower_global.z, "", "")
      }
      else if (string_cmp(t, "upper_global.z") || string_cmp(t, "zmax"))
      {
        GET_OR_CHOOSE_A_REAL(upper_global.z, "", "")
      }
      else if (string_cmp(t, "boundary_condition") || string_cmp(t, "bc"))
      {
        GET_OR_CHOOSE_A_INT_3D_VECTOR(boundary_condition, "", "")
      }
      else if (string_cmp(t, "me"))
      {
        std::cout << "ME: " << me << std::endl;
      }
      else if (string_cmp(t, "info"))
      {
        std::cout << "MPI Rank: " << me 
                  << " local.x [" << lower_local.x << " , " <<  upper_local.x << "]" 
                  << " local.y [" << lower_local.y << " , " <<  upper_local.y << "]" 
                  << " local.z [" << lower_local.z << " , " <<  upper_local.z << "]" 
                  << " global.x [" << lower_global.x << " , " <<  upper_global.x << "]" 
                  << " global.y [" << lower_global.y << " , " <<  upper_global.y << "]" 
                  << " global.z [" << lower_global.z << " , " <<  upper_global.z << "]" 
                  << std::endl;
      }
      else if (string_cmp(t, "generate"))
      {
        generate();
      }
      else
        error->all(FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
    }
    return in_file;
  }

  void Box::generate()
  {
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
    calculate_local_domain();
#else
    calculate_procs_grid();
    calculate_local_domain();
    std::string s = "boundary condition: ";
    s += std::to_string(boundary_condition.x) + " " + std::to_string(boundary_condition.y) + " " + std::to_string(boundary_condition.z);
    output->info(s);
#endif

    // it is defined for one-domain case.
    half_edge = 0.5 * (upper_global - lower_global); // XXX
  }

  void Box::calculate_local_domain()
  {
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
    lower_local.x = lower_global.x;
    lower_local.y = lower_global.y;
    lower_local.z = lower_global.z;

    upper_local.x = upper_global.x;
    upper_local.y = upper_global.y;
    upper_local.z = upper_global.z;
#else
    lower_local.x = lower_global.x + (upper_global.x - lower_global.x) * grid_index_x / nprocs_x;
    lower_local.y = lower_global.y + (upper_global.y - lower_global.y) * grid_index_y / nprocs_y;
    lower_local.z = lower_global.z + (upper_global.z - lower_global.z) * grid_index_z / nprocs_z;

    upper_local.x = lower_global.x + (upper_global.x - lower_global.x) * (grid_index_x + 1) / nprocs_x;
    upper_local.y = lower_global.y + (upper_global.y - lower_global.y) * (grid_index_y + 1) / nprocs_y;
    upper_local.z = lower_global.z + (upper_global.z - lower_global.z) * (grid_index_z + 1) / nprocs_z;
#endif
  }

  double Box::fix_distance_x(double d)
  {
#ifdef CAVIAR_WITH_MPI
    error->all(FC_FILE_LINE_FUNC, "not implemented.");
    return 0.0;
#else
    if (boundary_condition.x == 1)
    {
      if (d > +half_edge.x)
      {
        d -= half_edge.x * 2.0;
      }
      if (d < -half_edge.x)
      {
        d += half_edge.x * 2.0;
      }
    }
    return d;
#endif
  }

  double Box::fix_distance_y(double d)
  {
#ifdef CAVIAR_WITH_MPI
    error->all(FC_FILE_LINE_FUNC, "not implemented.");
    return 0.0;
#else
    if (boundary_condition.y == 1)
    {
      if (d > +half_edge.y)
      {
        d -= half_edge.y * 2.0;
      }
      if (d < -half_edge.y)
      {
        d += half_edge.y * 2.0;
      }
    }
    return d;
#endif
  }

  double Box::fix_distance_z(double d)
  {
#ifdef CAVIAR_WITH_MPI
    error->all(FC_FILE_LINE_FUNC, "not implemented.");
    return 0.0;
#else
    if (boundary_condition.z == 1)
    {
      if (d > +half_edge.z)
      {
        d -= half_edge.z * 2.0;
      }
      if (d < -half_edge.z)
      {
        d += half_edge.z * 2.0;
      }
    }
    return d;
#endif
  }

  caviar::Vector<double> Box::fix_distance(caviar::Vector<double> v)
  {
    return caviar::Vector<double>{fix_distance_x(v.x),
                                  fix_distance_y(v.y),
                                  fix_distance_z(v.z)};
  }

} // domain

CAVIAR_NAMESPACE_CLOSE
