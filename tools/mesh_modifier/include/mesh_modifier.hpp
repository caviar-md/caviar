
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

#ifndef MESHMODIFIER_H
#define MESHMODIFIER_H

#include <string>
#include <fstream>
#include <vector>

#include "unv_container.h"
#include "vector.h"
#include "point_condition.h"

namespace mesh_modifier
{

  class Mesh_modifier
  {

  public:
    Mesh_modifier();
    ~Mesh_modifier();

    void cat_by_string(const std::string &);
    void cat_by_line(const std::string &);
    void cat_by_line_into_token(const std::string &);

  public:
    void cat(const std::string &, bool unsupported);

  private:
    void cat_udn_ignore(std::ifstream &, int);
    void cat_udn_unsupported(std::ifstream &, int);
    void cat_udn_2411(std::ifstream &); // vertex
    void cat_udn_2412(std::ifstream &); // shape
    void cat_udn_2467(std::ifstream &); // group

  public:
    void import(const std::string &, bool unsupported);
    void import_vtk(const std::string &, bool unsupported);

  private:
    void import_udn_ignore(std::ifstream &, int);
    void import_udn_unsupported(std::ifstream &, int);
    void import_udn_2411(std::ifstream &); // vertex
    void import_udn_2412(std::ifstream &); // shape
    void import_udn_2467(std::ifstream &); // group

  public:
    void export_stl_boundary(const std::string &);

  private:
    void export_stl_headers(std::ofstream &);
    void export_stl_polygons(std::ofstream &);
    void export_stl_enders(std::ofstream &);

  public:
    void export_vtk(const std::string &);
    void export_vtk_boundary(const std::string &);

  private:
    void export_vtk_polygons(std::ofstream &);
    void export_vtk_headers(std::ofstream &);
    void export_vtk_points(std::ofstream &);
    void export_vtk_cells(std::ofstream &);
    void export_vtk_cell_types(std::ofstream &);
    void export_vtk_cell_data(std::ofstream &);

  public:
    void export_(const std::string &, bool unsupported);

  private:
    void export_udn_ignore(std::ofstream &);
    void export_udn_unsupported(std::ofstream &);
    void export_udn_2411(std::ofstream &); // vertex
    void export_udn_2412(std::ofstream &); // shape
    void export_udn_2467(std::ofstream &); // group

  public:
    void merge_all(bool unsupported);

  private:
    void add_to_unv_container(Unv_container &c1, Unv_container &c2);
    void fix_udn_2411_element_labels();
    void fix_udn_2412_element_labels();
    void fix_udn_2412_node_labels();
    void fix_udn_2467_element_labels();
    void fix_udn_2467_entity_tags();

  public:
    void remove_similar_points(const double tol_sqr);
    void remove_unnecessary_points(); // TODO
  private:
  public:
    void remove_hexa_internals(const double tol_sqr);
    unsigned count_quads();
    unsigned count_edges();

  private:
    void remove_hexa_internal_groups(std::vector<unsigned> &);
    void find_hexahedrons(std::vector<std::vector<unsigned>> &,
                          std::vector<std::vector<unsigned>> &);

  public:
    void convert_tetra_to_hexa();

  private:
    void find_all_tetra(std::vector<std::vector<unsigned>> &);
    unsigned index_of_point(unsigned label);

    unsigned add_point_to_2411(Vector<double> &);
    unsigned add_edge_to_2412(unsigned, unsigned);
    unsigned add_quad_to_2412(unsigned, unsigned, unsigned, unsigned);
    unsigned add_hexahedron_to_2412(unsigned, unsigned, unsigned, unsigned,
                                    unsigned, unsigned, unsigned, unsigned);
    void remove_edge_from_2412(unsigned, unsigned);
    void remove_triangle_from_2412(unsigned, unsigned, unsigned);
    void remove_all_tetra_edges(std::vector<unsigned> &);
    void remove_all_triangles();
    void remove_all_tetrahedrons();

  public:
    void merge_close_vertices(const double); // TODO
    void refine_mesh();                      // TODO
    void coarse_mesh();                      // TODO

    // XXX UTILITY XXX
  public:
    void add_face_to_group_with_condition(const std::string gr, const std::vector<Point_condition> *);
    void scale_vertices(const std::vector<Point_condition> *, const Vector<double> center, const Vector<double> sc);

  private:
    void remove_edge(std::vector<unsigned> &, std::vector<unsigned> &);
    void remove_quad(std::vector<unsigned> &, std::vector<unsigned> &);
    void add_edge_to_vector(std::vector<std::vector<unsigned>> &, unsigned, unsigned);
    void add_triangle_to_vector(std::vector<std::vector<unsigned>> &,
                                unsigned, unsigned, unsigned);
    void add_quad_to_vector(std::vector<std::vector<unsigned>> &,
                            unsigned, unsigned, unsigned, unsigned);

    void add_tetrahedron_to_vector(std::vector<std::vector<unsigned>> &,
                                   unsigned, unsigned, unsigned, unsigned);

    void extract_udn_2467_ids(std::vector<unsigned> &ids);
    void extract_udn_2467_index_ids(std::vector<unsigned> &ids);

    void remove_label_from_udn_2467(unsigned element_label,
                                    unsigned group_id);

    void add_label_to_udn_2467(unsigned element_label,
                               unsigned group_id);

    void make_hexa_edge_from_tetra_edge(std::vector<unsigned> &udn_2467_ids,
                                        std::vector<unsigned> &tetra_edge_index);
    void make_quad_from_triangle(std::vector<unsigned> &udn_2467_ids);
    void make_hexa_from_tetra(std::vector<unsigned> &udn_2467_ids);

  private:
    void make_label_to_index(const std::vector<Universal_dataset_number_2411> &,
                             std::vector<unsigned> &);
    void make_label_to_index(const std::vector<Universal_dataset_number_2412> &,
                             std::vector<unsigned> &);

  private:
    std::vector<Unv_container> unv_container;
  };

}

#endif
