
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

#include "caviar/objects/shape/polyhedron/format_vtk_reader.h"
#include "caviar/objects/shape/polyhedron/polyhedron.h"
#include "caviar/objects/shape/polyhedron/preprocess.h"
#include "caviar/utility/interpreter_io_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace shape
{
  namespace polyhedron
  {

    Format_vtk_reader::Format_vtk_reader(CAVIAR *fptr) : Pointers{fptr} {}

    Format_vtk_reader::~Format_vtk_reader() {}

    void Format_vtk_reader::read_polyhedron(shape::polyhedron::Polyhedron &p_object, const std::string &file)
    {
      class caviar::interpreter::Parser *parser = new caviar::interpreter::Parser{fptr, file};
      std::cout << "info: Vtk_file: reading a Vtk_file vtk file: " << file << std::endl;

      auto &vertex = p_object.vertex;
      auto &face = p_object.face;
      auto &edges = p_object.edges;
      auto &vertex_map = p_object.vertex_map;
      auto &face_id = p_object.face_id;

      if (face.size() != face_id.size())
      {
        error->all(FC_FILE_LINE_FUNC, "face.size() != face_id.size()");
      }

      std::cout << "info: Vtk_file: recording vtk file as a 3D shape with index "
                << " ? " << std::endl;

      while (true)
      {

        auto t0 = parser->get_val_token();
        struct caviar::interpreter::Token t1;
        if (t0.kind == caviar::interpreter::Kind::eof)
          break;

        if (t0.kind == caviar::interpreter::Kind::identifier)
        {
          if (t0.string_value == "POINTS")
          {
            output->info("Vtk_file: read vtk: POINTS");
            parser->get_literal_int(); // Number of points // useless
            parser->get_val_token();
            parser->end_of_line();
            int xyz = 0;
            Real_t x = 0, y = 0, z = 0;
            while (true)
            {
              t1 = parser->get_val_token();
              if (t1.kind == caviar::interpreter::Kind::identifier)
                break;
              if (t1.kind == caviar::interpreter::Kind::eol)
                t1 = parser->get_val_token();
              if (t1.kind == caviar::interpreter::Kind::identifier)
                break;
              if (xyz == 0)
              {
                if (t1.kind == caviar::interpreter::Kind::int_number)
                  x = t1.int_value;
                else
                  x = t1.real_value;
              }
              else if (xyz == 1)
              {
                if (t1.kind == caviar::interpreter::Kind::int_number)
                  y = t1.int_value;
                else
                  y = t1.real_value;
              }
              else
              {
                if (t1.kind == caviar::interpreter::Kind::int_number)
                  z = t1.int_value;
                else
                  z = t1.real_value;
                vertex.push_back(Vector<Real_t>{x, y, z});
              }
              ++xyz;
              if (xyz == 3)
                xyz = 0;
            }

            polyhedron::Preprocess p_pre(fptr);
            p_pre.merge_vertices(p_object);
          }
        }

        if (t1.string_value == "VERTICES")
        {
          output->info("Vtk_file: read vtk: VERTICES");
          while (true)
          {
            t1 = parser->get_val_token();
            if (t1.kind == caviar::interpreter::Kind::identifier)
            {
              // polyhedron::Preprocess p_pre (fptr);
              // p_pre.merge_vertices(p_object);
              break;
            }
          }
        }

        if (t1.string_value == "LINES")
        {
          output->info("Vtk_file: read vtk: LINES");
          while (true)
          {
            t1 = parser->get_val_token();
            if (t1.kind == caviar::interpreter::Kind::identifier)
            {
              //          std::cout<<"identifier " << t1.string_value<<std::endl;
              break;
            }
          }
        }

        if (t1.string_value == "POLYGONS" || t1.string_value == "CELLS")
        {
          //      faces_of_vertex.resize (vertex.size());
          output->info("Vtk_file: read vtk: POLYGONS");
          unsigned int no_polygons = parser->get_literal_int();
          parser->get_literal_int(); // useless variable
          parser->end_of_line();
          for (unsigned int i = 0; i < no_polygons; ++i)
          {
            parser->get_literal_int(); // useless variable
            unsigned int num1 = parser->get_literal_int();
            unsigned int num2 = parser->get_literal_int();
            unsigned int num3 = parser->get_literal_int();
            //        std::cout << num << " " << num1 << " " << num2 << " " << num3 << "\n";

            num1 = (vertex_map[num1].size() == 0) ? num1 : vertex_map[num1][1]; // uses vertex_map instead of vertex
            num2 = (vertex_map[num2].size() == 0) ? num2 : vertex_map[num2][1]; //
            num3 = (vertex_map[num3].size() == 0) ? num3 : vertex_map[num3][1]; //

            std::vector<unsigned int> gons;
            gons.push_back(num1);
            gons.push_back(num2);
            gons.push_back(num3);
            face.push_back(gons);
            face_id.push_back(-1); // invalid face_id;

            //        faces_of_vertex[num1].push_back(i); faces_of_vertex[num2].push_back(i); faces_of_vertex[num3].push_back(i);

            std::map<std::vector<unsigned int>, std::vector<unsigned int>>::iterator it_edges;
#define ADD_EDGE(NUM1, NUM2)                          \
  {                                                   \
    std::vector<unsigned int> check_face = {i};       \
    std::vector<unsigned int> temp_edge;              \
    if (NUM1 < NUM2)                                  \
      temp_edge = {NUM1, NUM2};                       \
    else                                              \
      temp_edge = {NUM2, NUM1};                       \
    it_edges = edges.find(temp_edge);                 \
    if (it_edges == edges.end())                      \
    {                                                 \
      edges.insert(make_pair(temp_edge, check_face)); \
    }                                                 \
    else                                              \
    {                                                 \
      it_edges->second.push_back(i);                  \
    }                                                 \
  }
            ADD_EDGE(num1, num2);
            ADD_EDGE(num2, num3);
            ADD_EDGE(num3, num1);
#undef ADD_EDGE
            parser->end_of_line();
          }

          break;
        }
      }

      delete parser;
    }

    void Format_vtk_reader::write_unstructured_vtk4(shape::polyhedron::Polyhedron &p_object,
                                                    const std::string st_out)
    {
      auto &vertex = p_object.vertex;
      auto &face = p_object.face;

      std::ofstream ofs(st_out.c_str());

      ofs << "# vtk DataFile Version 4.0\n"
             "vtk output\n"
             "ASCII\n"
             "DATASET UNSTRUCTURED_GRID\n"
             "POINTS "
          << vertex.size() << " float\n";

      for (auto &&v : vertex)
        ofs << v.x << " " << v.y << " " << v.z << "\n";

      ofs << "\n";

      int sum_face_vertex = 0;
      for (auto &&f : face)
        sum_face_vertex += f.size();

      ofs << "CELLS " << face.size() << " " << face.size() + sum_face_vertex << "\n";
      for (auto &&f : face)
      {
        ofs << f.size();
        for (auto &&i : f)
        {
          ofs << " " << i;
        }
        ofs << "\n";
      }
    }

    void Format_vtk_reader::write_polydata_vtk4(shape::polyhedron::Polyhedron &p_object,
                                                const std::string st_out)
    {
      auto &vertex = p_object.vertex;
      auto &face = p_object.face;

      std::ofstream ofs(st_out.c_str());

      ofs << "# vtk DataFile Version 4.0\n"
             "vtk output\n"
             "ASCII\n"
             "DATASET POLYDATA\n"
             "POINTS "
          << vertex.size() << " float\n";

      for (auto &&v : vertex)
        ofs << v.x << " " << v.y << " " << v.z << "\n";

      ofs << "\n";

      int sum_face_vertex = 0;
      for (auto &&f : face)
        sum_face_vertex += f.size();

      ofs << "POLYGONS " << face.size() << " " << face.size() + sum_face_vertex << "\n";
      for (auto &&f : face)
      {
        ofs << f.size();
        for (auto &&i : f)
        {
          ofs << " " << i;
        }
        ofs << "\n";
      }
    }

  } // polyhedron
} // shape

CAVIAR_NAMESPACE_CLOSE
