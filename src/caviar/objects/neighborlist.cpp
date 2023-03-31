
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

#include "caviar/objects/neighborlist.h"
#include "caviar/interpreter/error.h"

CAVIAR_NAMESPACE_OPEN

Neighborlist::Neighborlist (CAVIAR *fptr) : Pointers{fptr}, 
    atom_data{nullptr} {
  FC_OBJECT_INITIALIZE
}

Neighborlist::~Neighborlist () {}

void Neighborlist::verify_settings () {
  
}

Vector<int> Neighborlist::binlist_index (const Vector<double> &a) {
  error->all("NOT implemented yet");
  return Vector<int> {static_cast<int>(a.x),0,0};
}
int Neighborlist::neigh_bin_index (const Vector<double> &a) {
  error->all("NOT implemented yet");
  return a.x;
}


CAVIAR_NAMESPACE_CLOSE

