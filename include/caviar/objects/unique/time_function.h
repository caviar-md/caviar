
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

#ifndef CAVIAR_OBJECTS_UNIQUE_TIMEFUNCTION_H
#define CAVIAR_OBJECTS_UNIQUE_TIMEFUNCTION_H

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
   * This class defines a function of time that can be used in other objects such as force_fields
   *
   */
  class Time_function : public Unique
  {
  public:
    Time_function();
    Time_function(class CAVIAR *);
    ~Time_function();
    bool read(caviar::interpreter::Parser *);
    void generate_export_file();
    void generate_formula();
    void verify_settings();
    double value() { return current_value; };
    void update_time_variable(double t);
    void calculate();

    std::string function_definition;
    double time_variable;
    double current_value;
    bool export_values_to_file;
    bool export_file_append;
    std::string export_file_name;
    std::ofstream ofs_time_value;

    mu::Parser *muParser;
  };

} // unique

CAVIAR_NAMESPACE_CLOSE

#endif
