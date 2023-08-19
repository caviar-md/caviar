
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// Implementation of Fix_bond is done by Nasrin Eyvazi
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

#include "caviar/objects/force_field/fix_bond.h"
#include "caviar/objects/atom_data/utility/bond.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Fix_bond::Fix_bond(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    type_i = -1;
    type_j = -1;
    btype = -1;
    blength = -1;
    Rmin = -1;
  }

  bool Fix_bond::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "type_i"))
      {
        GET_OR_CHOOSE_A_INT(type_i, "", "")
        if (type_i < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "type_i have to be non-negative integer.");
      }
      else if (string_cmp(t, "type_j"))
      {
        GET_OR_CHOOSE_A_INT(type_j, "", "")
        if (type_j < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "type_j have to be non-negative integer.");
      }
      else if (string_cmp(t, "btype"))
      {
        GET_OR_CHOOSE_A_INT(btype, "", "")
        if (btype < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "bond type have to be non-negative integer.");
      }
      else if (string_cmp(t, "blength"))
      {
        GET_OR_CHOOSE_A_REAL(blength, "", "")
        if (blength < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "bond length have to non-negative.");
      }
      else if (string_cmp(t, "Rmin"))
      {
        GET_OR_CHOOSE_A_REAL(Rmin, "", "")
        if (Rmin < 0.0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "R min have to non-negative.");
      }
      else if (string_cmp(t, "bond_limit"))
      {
        GET_A_STDVECTOR_INT_ELEMENT(bond_limit)
        if (vector_value < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Bond limit have to be non-negative.");
      }
      else if (string_cmp(t, "set_neighborlist") || string_cmp(t, "neighborlist"))
      {
        FIND_OBJECT_BY_NAME(neighborlist, it)
        neighborlist = object_container->neighborlist[it->second.index];
      }
      else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Fix_bond::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    if (type_i == -1)
      error->all(FC_FILE_LINE_FUNC, "(type_i == -1)");
    if (type_j == -1)
      error->all(FC_FILE_LINE_FUNC, "(type_j == -1)");
    if (btype == -1)
      error->all(FC_FILE_LINE_FUNC, "(btype == -1)");
    if (blength < 0)
      error->all(FC_FILE_LINE_FUNC, "(blength < 0)");
    if (Rmin < 0)
      error->all(FC_FILE_LINE_FUNC, "(Rmin < 0)");

    auto type_max = type_i > type_j ? type_i : type_j;
    if ((int)bond_limit.size() < type_max + 1)
      bond_limit.resize(type_max + 1, 0);
    my_mpi_rank = atom_data->get_mpi_rank();
  }

  void Fix_bond::create_atomic_bond()
  {

    auto &pos = atom_data->atom_struct_owned.position;
    auto &bond_count = atom_data->atom_struct_owned.atomic_bond_count;

    for (unsigned int i = 0; i < pos.size(); ++i)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      for (unsigned int j = i + 1; j < pos.size(); ++j) // If the other particle is ghost, new bond will not be created.
      {
#ifdef CAVIAR_WITH_MPI
        if (atom_data->atom_struct_owned.mpi_rank[j] != my_mpi_rank) // If the other particle is ghost, new bond will not be created.
          continue;
#endif
        if ((int)atom_data->atom_struct_owned.type[i] == type_i && (int)atom_data->atom_struct_owned.type[j] == type_j)
        {
          ;
          const auto dr = pos[j] - pos[i];
          const auto dr_sq = dr * dr;

          if (dr_sq < Rmin * Rmin)
          {

            atom_data::Bond b;
            b.type = btype;
            b.length = blength;
            b.id_1 = i;
            b.id_2 = j;

            if (bond_count[i] < bond_limit[type_i] && bond_count[j] < bond_limit[type_j])
            {
              if (!atom_data->check_atomic_bond_exist(b))
              {
                atom_data->add_atomic_bond(b);
              }
            }
          }
        }
      }
    }
  }
  void Fix_bond::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    create_atomic_bond();
  }
} // force_field

CAVIAR_NAMESPACE_CLOSE
