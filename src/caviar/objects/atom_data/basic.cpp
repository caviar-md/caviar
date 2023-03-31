
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

#include "caviar/objects/atom_data/basic.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/unique/atom.h"
#include "caviar/objects/unique/atom_group.h"
#include "caviar/objects/unique/atom_list.h"
#include "caviar/objects/unique/molecule.h"
#include "caviar/objects/unique/molecule_group.h"
#include "caviar/objects/unique/molecule_list.h"
#include "caviar/objects/neighborlist/cell_list.h"
#include "caviar/interpreter/object_handler/preprocessors_new.h"

namespace caviar {

namespace atom_data {

Basic::Basic (CAVIAR *fptr) : Atom_data{fptr} {
  FC_OBJECT_INITIALIZE_INFO
}

Basic::~Basic() {}

} //atom_data

} // namespace caviar

