
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

#include "caviar/objects/force_field.hpp"

namespace caviar
{

  namespace force_field
  {

    /**
     * This class has electrostatic force-field.
     *
     *
     */
    class Electrostatic : public Force_field
    {
    public:
      Electrostatic(class CAVIAR *);
      ~Electrostatic() {};
      double potential(const Vector3d<double> &);
      double potential(const int);

      Vector3d<double> field(const Vector3d<double> &);
      Vector3d<double> field(const int);

      double energy();

      bool read(class caviar::interpreter::Parser *);
      void verify_settings();
      void calculate_acceleration();

    public:
      // std::vector<std::vector<double>> epsilon, sigma;
      std::vector<std::vector<double>> lambda;
      bool lambda_is_set = false;
      double k_electrostatic;
      Vector3d<double> external_field;
    };

  } // force_field

}
