
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

#include "caviar/objects/constraint.hpp"
#include "caviar/CAVIAR.hpp"

namespace caviar
{

  Constraint::Constraint(CAVIAR *fptr) : caviar_{fptr},
                                   comm{fptr->comm},
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
                                   err_flag{fptr->err_flag} {
                                   atom_data=nullptr;
                                                             FC_OBJECT_INITIALIZE}

                                         Constraint::~Constraint()
  {
  }

  void Constraint::verify_settings()
  {
  }

  void Constraint::apply(int64_t) {}
  void Constraint::apply_shake(int64_t) {}
  void Constraint::fix_position(int64_t) {}
  void Constraint::fix_velocity(int64_t, bool &) {}
  void Constraint::apply_barostat(int64_t, bool &) {}
  void Constraint::apply_thermostat(int64_t, bool &) {}
  void Constraint::fix_acceleration(int64_t) {}

}
