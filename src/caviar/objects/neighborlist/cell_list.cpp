
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

#include "caviar/objects/neighborlist/cell_list.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/domain.h"

#include <cmath>

namespace caviar {
namespace objects {
namespace neighborlist {
inline int int_floor(double x) 
{ 
    return (int)(x+100000) - 100000; 
}

Cell_list::Cell_list (CAVIAR *fptr) : Neighborlist{fptr}, domain{nullptr}
{
  FC_OBJECT_INITIALIZE_INFO
  cutoff = 0;
  cutoff_neighlist = 0;
  make_neighlist = false;
}

bool Cell_list::read (caviar::interpreter::Parser *parser) {
  FC_OBJECT_READ_INFO
  bool in_file = true;
  while(true) {
    GET_A_TOKEN_FOR_CREATION
    auto t = token.string_value;
    if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
      FIND_OBJECT_BY_NAME(atom_data,it)
      atom_data = object_container->atom_data[it->second.index];
    } else if (string_cmp(t,"cutoff")) {
      GET_OR_CHOOSE_A_REAL(cutoff,"","")
      if (cutoff < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "cutoff have to non-negative."); 
    } else if (string_cmp(t,"cutoff_neighlist")) {
      GET_OR_CHOOSE_A_REAL(cutoff_neighlist,"","")
      if (cutoff_neighlist < 0.0) error->all (FC_FILE_LINE_FUNC_PARSE, "cutoff_neighlist have to non-negative."); 
    } else if (string_cmp(t,"make_neighlist")) {
      make_neighlist = true;
    } else if (string_cmp(t,"set_domain") || string_cmp(t,"domain")) {
      FIND_OBJECT_BY_NAME(domain,it)
      domain = object_container->domain[it->second.index];
    } else error->all (FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
  }
  return in_file;
}

void Cell_list::init () {

  FC_NULLPTR_CHECK(atom_data)
  FC_NULLPTR_CHECK(domain)
  if (make_neighlist == false)
    output->warning("'make_neighlist' is false. You can only use force_field "
                    "classes that have 'cell_list' implementation. if not, you "
                    "may get a 'segmentation fault' error.");

  auto dd = domain->upper_local - domain->lower_local; // NO DOMAIN DECOMPOSITION
  auto bc = domain->boundary_condition;


  no_bins.x = std::ceil(dd.x/cutoff) ;
  no_bins.y = std::ceil(dd.y/cutoff) ;
  no_bins.z = std::ceil(dd.z/cutoff) ;

  if (bc.x==1) no_bins.x  += 2;
  if (bc.y==1) no_bins.y  += 2;
  if (bc.z==1) no_bins.z  += 2;

  std::string s = "no_bins.x : " + std::to_string(no_bins.x)
                + " ,no_bins.y : " + std::to_string(no_bins.y)
                + " ,no_bins.z : " + std::to_string(no_bins.z);
  output->info (s);

  binlist.resize(no_bins.x);

  for (unsigned i=0 ; i < binlist.size(); ++i)
    binlist[i].resize(no_bins.y);

  for (unsigned i=0 ; i < binlist.size(); ++i)
    for (unsigned j=0 ; j < binlist[i].size(); ++j)
      binlist[i][j].resize(no_bins.z);

  make_neigh_bin();
}

bool Cell_list::rebuild_neighlist () {
  return true;
}

Vector<int> Cell_list::binlist_index (const Vector<double> & pos) {
  auto ind = domain->boundary_condition;
  ind.x += int_floor( (pos.x - domain->lower_local.x)/cutoff );
  ind.y += int_floor( (pos.y - domain->lower_local.y)/cutoff );
  ind.z += int_floor( (pos.z - domain->lower_local.z)/cutoff );
  if (ind.x < 0) ind.x = 0;
  if (ind.y < 0) ind.y = 0;
  if (ind.z < 0) ind.z = 0;
  if (ind.x >= no_bins.x) ind.x = no_bins.x - 1;
  if (ind.y >= no_bins.y) ind.y = no_bins.y - 1;
  if (ind.z >= no_bins.z) ind.z = no_bins.z - 1;
  return ind;
}

int Cell_list::neigh_bin_index (const Vector<double> & pos) {
  const auto ind = binlist_index (pos);
  return ind.x + no_bins.x*ind.y + no_bins.x*no_bins.y*ind.z;
}

void Cell_list::build_binlist () {
  for (auto && bx : binlist)
    for (auto && by : bx)
      for (auto && bz : by)
        bz.clear();

  const auto &pos = atom_data->owned.position;
  const auto &pos_ghost = atom_data->ghost.position;
  const auto pos_size = pos.size();

  for (unsigned int i=0; i<pos_size; ++i) {
    const auto ind = binlist_index (pos[i]);
    binlist[ind.x][ind.y][ind.z].emplace_back(i);
  }

  for (unsigned int i=0; i<pos_ghost.size(); ++i) {
    const auto ind = binlist_index (pos_ghost[i]);
    binlist[ind.x][ind.y][ind.z].emplace_back(i+pos_size);
  }
}


void Cell_list::build_neighlist () {
  build_binlist ();
  if (make_neighlist) {
    const auto &pos = atom_data->owned.position;
    const int pos_size = pos.size();
    const auto cutoff_neighlist_sq = cutoff_neighlist*cutoff_neighlist;
    const auto &nb = neigh_bin;
    neighlist.clear ();
    neighlist.resize (pos_size);


    for (int i=0; i<pos_size; ++i) {

      auto nb_i = neigh_bin_index (pos[i]);

      for (unsigned nb_j = 0; nb_j < nb[nb_i].size(); ++nb_j) {
        const auto &nb_ij = nb[nb_i][nb_j];
        for (unsigned j = 0; j < binlist [nb_ij.x] [nb_ij.y] [nb_ij.z].size(); ++j) {

          const auto ind_j = binlist[nb_ij.x] [nb_ij.y] [nb_ij.z][j];
          if (ind_j>i) {
            auto dr = pos[i];
            if (ind_j >= pos_size) // Ghost condition
              dr -= atom_data->ghost.position[ind_j-pos_size];
            else 
              dr -= atom_data->owned.position[ind_j]; 
            if (dr*dr < cutoff_neighlist_sq) 
              neighlist[i].emplace_back (ind_j);
          }
        }
      }
    }
  }
}


// this function will be called only once at the start or after any change in
// the Cell_list.
// Here we make the 'neigh_bin'. It contains the list of neighbors (of cell type)
// any cell
// could have, including itself.
void Cell_list::make_neigh_bin () {
  neigh_bin.clear();
  neigh_bin.resize(no_bins.x*no_bins.y*no_bins.z);

  Vector<int> ind {0,0,0};

  for (ind.x = 0; ind.x < no_bins.x; ++ind.x) {
    const int ind_x_min = ind.x > 0 ? ind.x - 1 : 0 ;
    const int ind_x_max = ind.x < no_bins.x - 1 ? ind.x + 1 : no_bins.x - 1 ;

    for (ind.y = 0; ind.y < no_bins.y; ++ind.y) {
      const int ind_y_min = ind.y > 0 ? ind.y - 1 : 0 ;
      const int ind_y_max = ind.y < no_bins.y - 1 ? ind.y + 1 : no_bins.y - 1 ;

      for (ind.z = 0; ind.z < no_bins.z; ++ind.z) {
        const int ind_z_min = ind.z > 0 ? ind.z - 1 : 0 ;
        const int ind_z_max = ind.z < no_bins.z - 1 ? ind.z + 1 : no_bins.z - 1 ;

        const int k = ind.x + no_bins.x*ind.y + no_bins.x*no_bins.y*ind.z;

        for (int ind_x_i = ind_x_min; ind_x_i <= ind_x_max; ++ind_x_i)
        for (int ind_y_i = ind_y_min; ind_y_i <= ind_y_max; ++ind_y_i)
        for (int ind_z_i = ind_z_min; ind_z_i <= ind_z_max; ++ind_z_i)      
          neigh_bin[k].emplace_back(ind_x_i, ind_y_i, ind_z_i);
        
      }
    }
  }
}

} //neighborlist
} //objects
} // namespace caviar

