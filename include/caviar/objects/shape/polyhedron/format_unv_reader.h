
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

#ifndef CAVIAR_OBJECTS_SHAPE_POLYHEDRON_FORMAT_UNV_READER_h
#define CAVIAR_OBJECTS_SHAPE_POLYHEDRON_FORMAT_UNV_READER_h

#include "caviar/utility/pointers.h"

CAVIAR_NAMESPACE_OPEN

namespace shape
{
  namespace polyhedron
  {
    struct Polyhedron;
    class Format_unv_reader : public Pointers
    {
    public:
      Format_unv_reader(class CAVIAR *);
      ~Format_unv_reader();

      void read_polyhedron(shape::polyhedron::Polyhedron &, const std::string &);

      void import_udn_ignore(std::ifstream &ifs, int udn_code);
      void import_udn_2411(shape::polyhedron::Polyhedron &p_object, std::ifstream &); // vertex
      void import_udn_2412(shape::polyhedron::Polyhedron &p_object, std::ifstream &); // shape
      void import_udn_2467(shape::polyhedron::Polyhedron &p_object, std::ifstream &); // group

      void write_unv(shape::polyhedron::Polyhedron &, const std::string st_out = "o_shape.unv");
      void export_udn_2411(shape::polyhedron::Polyhedron &p_object, std::ofstream &ofs);
      void export_udn_2412(shape::polyhedron::Polyhedron &p_object, std::ofstream &ofs);
      void export_udn_2467(shape::polyhedron::Polyhedron &p_object, std::ofstream &ofs);

      void make_label_to_index(const std::vector<int> &u, std::vector<int> &v);

      std::vector<int> node_label;
      std::vector<int> node_label_to_index;
      std::vector<int> face_label;
      std::vector<int> face_label_to_index;
    };
  } // polyhedron
} // shape

CAVIAR_NAMESPACE_CLOSE

#endif
