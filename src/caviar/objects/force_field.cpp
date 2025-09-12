
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

#include "caviar/objects/force_field.hpp"
#include "caviar/interpreter/error.hpp"
#include "caviar/CAVIAR.hpp"

namespace caviar
{

  Force_field::Force_field(CAVIAR *fptr) : caviar_{fptr}, comm{fptr->comm},
                                           error{fptr->error},
                                           output{fptr->output},
                                           input{fptr->input},
                                           object_handler{fptr->object_handler},
                                           object_container{fptr->object_container},
                                           object_creator{fptr->object_creator},
                                           log{fptr->log},
                                           in{fptr->in},
                                           out{fptr->out},
                                           err{fptr->err},
                                           log_flag{fptr->log_flag},
                                           out_flag{fptr->out_flag},
                                           err_flag{fptr->err_flag}

  {
    FC_OBJECT_INITIALIZE
  }

  void Force_field::verify_settings()
  {
  }

  double Force_field::energy()
  {
    error->all(FC_FILE_LINE_FUNC, "The energy calculation of this force_field is not implemented");
    return 0.0;
  }

  double Force_field::potential(const Vector<double> &)
  {
    error->all(FC_FILE_LINE_FUNC, "The potential calculation of this force_field is not implemented");
    return 0.0;
  }

  double Force_field::potential(const int)
  {
    error->all(FC_FILE_LINE_FUNC, "The potential calculation of this force_field is not implemented");
    return 0.0;
  }

  Vector<double> Force_field::field(const Vector<double> &)
  {
    error->all(FC_FILE_LINE_FUNC, "The field calculation of this force_field is not implemented");
    return Vector<double>{0, 0, 0};
  }

  Vector<double> Force_field::field(const int)
  {
    error->all(FC_FILE_LINE_FUNC, "The field calculation of this force_field is not implemented");
    return Vector<double>{0, 0, 0};
  }

  void Force_field::scale_position(double, caviar::Vector<int>)
  {
    error->all(FC_FILE_LINE_FUNC, "The scale_position of this force_field is not implemented");
  }

}
