
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

#pragma once

#include "caviar/utility/objects_common_headers.hpp"

namespace caviar
{

  class Atom_data;
  class Domain;
  class Neighborlist;

  /**
   * This class is the base class for all the force-fields.
   * The  child object does not need to have all of the virtual function implemented.
   * Except the abstract ones.
   */
  class Force_field //
  {
  public:
    Force_field(class CAVIAR *);
    virtual ~Force_field() {};
    virtual bool read(class caviar::interpreter::Parser *) = 0;
    virtual void calculate_acceleration() = 0;
    virtual double energy();
    virtual double potential(const Vector3d<double> &);
    virtual double potential(const int);
    virtual Vector3d<double> field(const Vector3d<double> &);
    virtual Vector3d<double> field(const int);

    /**
     * Used in barostat scaling for geometrical forces.
     *
     */
    virtual void scale_position(double scale_ratio, caviar::Vector3d<int> scale_axis);

    double cutoff;
    class Atom_data *atom_data = nullptr;
    class Domain *domain = nullptr;
    class Neighborlist *neighborlist = nullptr;
    /**
     * MPI rank of the classs
     */
    int my_mpi_rank = -1;
    FC_BASE_OBJECT_COMMON_TOOLS

  };

}
