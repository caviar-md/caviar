
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


#ifndef CAVIAR_OBJECTS_FORCEFIELD_H
#define CAVIAR_OBJECTS_FORCEFIELD_H

#include "caviar/utility/objects_common_headers.h"

namespace caviar {
namespace objects {
class Atom_data; 
class Domain;
class Neighborlist; 

/**
 * This class is the base class for all the force-fields.
 * The  child object does not need to have all of the virtual function implemented.
 * Except the abstract ones.
 */
class Force_field  : public Pointers {
 public:
  Force_field (class CAVIAR *);
  virtual ~Force_field () {};
  virtual bool read (class caviar::interpreter::Parser *) = 0;
  virtual void calculate_acceleration () = 0;
  virtual double energy();
  virtual double potential (const Vector<double> &);
  virtual double potential (const int);
  virtual Vector<double> field (const Vector<double> &);
  virtual Vector<double> field (const int);
  double cutoff;
  class objects::Atom_data *atom_data;
  class objects::Domain *domain;
  class objects::Neighborlist *neighborlist;

  FC_BASE_OBJECT_COMMON_TOOLS
};

} //objects

} // namespace caviar

#endif
