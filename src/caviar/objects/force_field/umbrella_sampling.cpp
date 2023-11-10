
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

#include "caviar/objects/force_field/umbrella_sampling.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/utility/file_utility.h"
#include "caviar/objects/neighborlist.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/domain.h"

#include <cmath>
#include <iomanip>

CAVIAR_NAMESPACE_OPEN

namespace force_field
{

  Umbrella_sampling::Umbrella_sampling(CAVIAR *fptr) : Force_field{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    fix_file_prefix();
  }

  bool Umbrella_sampling::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "elastic_coef"))
      {
        GET_OR_CHOOSE_A_REAL(elastic_coef, "", "")
        if (elastic_coef < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "Elastic coef. have to be non-negative.");
      }
      else if (string_cmp(t, "atom_id"))
      {
        GET_OR_CHOOSE_A_INT(atom_id, "", "")
        if (atom_id < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "atom_id have to be non-negative.");
      }
        else if (string_cmp(t, "step"))
      {
        GET_OR_CHOOSE_A_INT(step, "", "")
        if (step <= 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "step have to be a positive number.");
      }
      else if (string_cmp(t, "set_position"))
      {
        double x = 0,y = 0,z = 0;
        GET_OR_CHOOSE_A_REAL(x, "", "")
        GET_OR_CHOOSE_A_REAL(y, "", "")
        GET_OR_CHOOSE_A_REAL(z, "", "")
        position = Vector<Real_t>{x, y, z};
      }
        else if (string_cmp(t, "set_position_x"))
      {
        GET_OR_CHOOSE_A_REAL(position.x, "", "")
      }
        else if (string_cmp(t, "set_position_y"))
      {
        GET_OR_CHOOSE_A_REAL(position.y, "", "")
      }
      
        else if (string_cmp(t, "set_position_z"))
      {
        GET_OR_CHOOSE_A_REAL(position.z, "", "")
      }
        else if (string_cmp(t, "set_temperature"))
      {
        GET_OR_CHOOSE_A_REAL(temperature, "", "")
        temperature_is_set = true;
      }
        else if (string_cmp(t, "reaction_coordinate"))
      {
        std::string s = "x";
        GET_A_STRING(s,"","");
        if (s == "x") reaction_coordinate = 'x';
        else if (s == "y") reaction_coordinate = 'y';
        else if (s == "z") reaction_coordinate = 'z';
        else
          error->all(FC_FILE_LINE_FUNC_PARSE,"Expected x, y or z as reaction coordinate");
      }
        else if (string_cmp(t, "file_prefix"))
      {
        auto t2 = parser->get_val_token();
        file_prefix = t2.string_value;
        fix_file_prefix();
      }
        else if (string_cmp(t, "init_production"))
      {
        init_production();
      }
        else if (string_cmp(t, "finish_production"))
      {
        finish_production();
      }
        else if (string_cmp(t, "init_metadata"))
      {
        init_metadata();
      }
        else if (string_cmp(t, "finish_metadata"))
      {
        finish_metadata();
      }
      else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else if (string_cmp(t, "set_domain") || string_cmp(t, "domain"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }

    return in_file;
  }

  void Umbrella_sampling::fix_file_prefix()
  {
    if (!is_absolute_path(file_prefix))
    {
      if (file_prefix == "")
        file_prefix = get_current_directory();
      else
        file_prefix = join_path(get_current_directory(), file_prefix);
    }
  }

  void Umbrella_sampling::init_metadata()
  {
    metadata_mode = true;

    std::string file_name_tmp = join_path(file_prefix, file_name_metadata + ".txt");
    ofs_metadata.open(file_name_tmp.c_str());

    std::string file_name_stat_tmp = join_path(file_prefix, file_name_stat + ".txt");
    ofs_stat.open(file_name_stat_tmp.c_str());
  }

  void Umbrella_sampling::finish_metadata()
  {
    metadata_mode = false;
    ofs_metadata.close();
    ofs_stat.close();
  }

  void Umbrella_sampling::init_production()
  {
    production_mode = true;
    step_counter = -1;
    std::string file_name_tmp = join_path(file_prefix, file_name_data + std::to_string(file_counter) + ".txt");
    ofs_data.open(file_name_tmp.c_str());
    if (metadata_mode == false)
    {
      error->all(FC_FILE_LINE_FUNC,"Expected metadata activation, i.e., 'init_metadata' call");
    }

    session_data.clear();
    position_mean_count = 0;
    position_mean = 0;

    ofs_metadata << file_name_tmp;
    position_min = 1e20;
    position_max = -1e20;
    if (reaction_coordinate == 'x')
    {
      ofs_metadata << " " << position.x;
    }
    else if (reaction_coordinate == 'y')
    {
      ofs_metadata << " " << position.y;
    }
    else if (reaction_coordinate == 'z')
    {
      ofs_metadata << " " << position.z;
    }
    else
    {
      error->all(FC_FILE_LINE_FUNC,"Expected x, y or z as reaction coordinate");
    }
    ofs_metadata << " " << elastic_coef ;
    ofs_metadata << " " << "0" ; // coecorrelation time for your time series used in monte-carlo

    if (temperature_is_set)
    {
      ofs_metadata << " " << temperature;
    }
    ofs_metadata << "\n";
  }

  void Umbrella_sampling::production_function()
  {
    if (!production_mode)
      return;

    step_counter++;

    if (step_counter % step != 0)
      return;

    int i = atom_data->atom_id_to_index[atom_id];
    auto p_i = atom_data->atom_struct_owned.position[i];

    ofs_data << step_counter;

    
    if (reaction_coordinate == 'x')
    {
      ofs_data << " " << p_i.x;
      if (position_min > p_i.x) position_min = p_i.x;
      if (position_max < p_i.x) position_max = p_i.x;
      position_mean += p_i.x;
    }
    else if (reaction_coordinate == 'y')
    {
      ofs_data << " " << p_i.y;
      if (position_min > p_i.y) position_min = p_i.y;
      if (position_max < p_i.y) position_max = p_i.y;
      position_mean += p_i.y;

    }
    else if (reaction_coordinate == 'z')
    {
      ofs_data << " " << p_i.z;
      if (position_min > p_i.z) position_min = p_i.z;
      if (position_max < p_i.z) position_max = p_i.z;
      position_mean += p_i.z;
    }
    else
    {
      error->all(FC_FILE_LINE_FUNC,"Expected x, y or z as reaction coordinate");
    }
    position_mean_count++;
    ofs_data << "\n";
    
    //
    //ofs_data << atom_data->atom_struct_owned.position[i] << " " << dr << " " << 0.5 * elastic_coef * dr * dr << "\n";
  }

  void Umbrella_sampling::finish_production()
  { 
    if (position_mean_count > 0)
      position_mean /= position_mean_count;
    
    ofs_stat << file_counter;
    if (reaction_coordinate == 'x')
    {
      ofs_stat << " " << position.x;
    }
    else if (reaction_coordinate == 'y')
    {
      ofs_stat << " " << position.y;
    }
    else if (reaction_coordinate == 'z')
    {
      ofs_stat << " " << position.z;
    }
    ofs_stat << " " << position_mean << " " << position_min << " " << position_max << "\n" << std::flush;

    ofs_data << std::flush;
    production_mode = false;
    ofs_data.close();
    file_counter++;

  }

  void Umbrella_sampling::verify_settings()
  {
    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(domain)
    //my_mpi_rank = atom_data->get_mpi_rank();
    if ((unsigned int) atom_id > atom_data->atom_id_to_index.size() - 1)
      error->all(FC_FILE_LINE_FUNC, "Invalid atom_id:" + std::to_string(atom_id));

    //int i = atom_data->atom_id_to_index[atom_id];

  }

  void Umbrella_sampling::calculate_acceleration()
  {
    FC_OBJECT_VERIFY_SETTINGS

    auto &type = atom_data->atom_struct_owned.type;
    auto &mass_inv = atom_data->atom_type_params.mass_inv;
    int i = atom_data->atom_id_to_index[atom_id];

    auto pos_i = atom_data->atom_struct_owned.position[atom_id];
    
#if defined(CAVIAR_WITH_MPI)
    dr = position - pos_i;
#else
    dr = domain->periodic_distance(position - pos_i);
#endif

    //const auto dr_sq = dr * dr;
    //const auto dr_norm = std::sqrt(dr_sq);
    //const auto dr_vec = dr / dr_norm;
    const auto force = -elastic_coef *dr ;

    atom_data->atom_struct_owned.acceleration[i] -= force * mass_inv[type[i]];

    production_function();

  }

} // force_field

CAVIAR_NAMESPACE_CLOSE
