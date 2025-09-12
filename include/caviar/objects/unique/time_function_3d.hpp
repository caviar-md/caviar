
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

#ifndef CAVIAR_OBJECTS_UNIQUE_TIMEFUNCTION3D_H
#define CAVIAR_OBJECTS_UNIQUE_TIMEFUNCTION3D_H

#include "caviar/objects/unique.h"
#include <random>

namespace mu
{
  class Parser;
}

CAVIAR_NAMESPACE_OPEN
class Parser;

namespace unique
{

  /**
   * This class defines a 3D function of time that can be used in other objects such as force_fields
   *
   */
  class Time_function_3d : public Unique
  {
  public:
    Time_function_3d();
    Time_function_3d(class CAVIAR *);
    ~Time_function_3d();
    bool read(caviar::interpreter::Parser *);
    void generate_export_file();
    void generate_formula();
    void verify_settings();
    Vector<double> value() { return current_value; };
    void update_time_variable(double t);
    void calculate();

    std::string function_definition_x = "0";
    std::string function_definition_y = "0";
    std::string function_definition_z = "0";
    double time_variable;
    Vector<double> current_value;
    bool export_values_to_file;
    bool export_file_append;
    std::string export_file_name;
    std::ofstream ofs_time_value;

    mu::Parser *muParser_x;
    mu::Parser *muParser_y;
    mu::Parser *muParser_z;
  };

} // unique

CAVIAR_NAMESPACE_CLOSE

#endif
