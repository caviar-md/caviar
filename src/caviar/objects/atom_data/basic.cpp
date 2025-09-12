
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

#include "caviar/objects/atom_data/basic.hpp"
#include "caviar/utility/interpreter_io_headers.hpp"
#include "caviar/objects/unique/atom.hpp"
#include "caviar/objects/unique/atom_group.hpp"
#include "caviar/objects/unique/atom_list.hpp"
#include "caviar/objects/unique/molecule.hpp"
#include "caviar/objects/unique/molecule_group.hpp"
#include "caviar/objects/unique/molecule_list.hpp"
#include "caviar/interpreter/object_handler/preprocessors_new.hpp"

CAVIAR_NAMESPACE_OPEN

namespace atom_data
{

  Basic::Basic(CAVIAR *fptr) : Atom_data{fptr} {
                                   FC_OBJECT_INITIALIZE_INFO}

                               Basic::~Basic()
  {
  }

} // atom_data

CAVIAR_NAMESPACE_CLOSE
