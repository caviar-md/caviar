
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

#include "caviar/objects/force_field/umbrella_sampling_g.hpp"
#include "caviar/utility/interpreter_io_headers.hpp"
#include "caviar/utility/file_utility.hpp"
#include "caviar/objects/neighborlist.hpp"
#include "caviar/objects/atom_data.hpp"
#include "caviar/objects/domain.hpp"
#include "caviar/objects/writer/atom_data.hpp"

#include <cmath>
#include <iomanip>

namespace caviar
{

  namespace force_field
  {

    Umbrella_sampling_g::Umbrella_sampling_g(CAVIAR *fptr) : Force_field{fptr}
    {
      FC_OBJECT_INITIALIZE_INFO
      fix_file_prefix();
      
    }

    Umbrella_sampling_g::~Umbrella_sampling_g() {
      if (writerXYZ != nullptr)
        delete writerXYZ;
      }

    bool Umbrella_sampling_g::read(caviar::interpreter::Parser *parser)
    {
      FC_OBJECT_READ_INFO
      bool in_file = true;

      while (true)
      {
        GET_A_TOKEN_FOR_CREATION
        auto t = token.string_value;
        FC_OBJECT_READ_INFO_STR
        if (string_cmp(t, "elastic_coef"))
        {
          GET_OR_CHOOSE_A_REAL(elastic_coef, "", "")
          if (elastic_coef < 0)
            error->all(FC_FILE_LINE_FUNC_PARSE, "Elastic coef. have to be non-negative.");
        }
        else if (string_cmp(t, "add_fixed_atom_id"))
        {
          int id = -1;
          GET_OR_CHOOSE_A_INT(id, "", "")
          if (id < 0)
            error->all(FC_FILE_LINE_FUNC_PARSE, "atom_id have to be non-negative.");

          if ((unsigned) id >= atom_data->atom_struct_owned.molecule_index.size())
            error->all(FC_FILE_LINE_FUNC_PARSE, "atom_id is larger than the number of atoms");

          fixed_atoms_id.push_back(id);
          fixed_atoms_resting_position.push_back(atom_data->atom_struct_owned.position[id]);

        }
        else if (string_cmp(t, "add_moving_atom_id"))
        {
          int id = -1;
          GET_OR_CHOOSE_A_INT(id, "", "")
          if (id < 0)
            error->all(FC_FILE_LINE_FUNC_PARSE, "atom_id have to be non-negative.");

          if ((unsigned) id >= atom_data->atom_struct_owned.molecule_index.size())
            error->all(FC_FILE_LINE_FUNC_PARSE, "atom_id is larger than the number of atoms");

          moving_atoms_id.push_back(id);
          moving_atoms_resting_position.push_back(atom_data->atom_struct_owned.position[id]);
          moving_atoms_initial_position.push_back(atom_data->atom_struct_owned.position[id]);

        }
          else if (string_cmp(t, "fixed_atom_id"))
        {
          GET_OR_CHOOSE_A_INT(fixed_atom_id, "", "")
          if (fixed_atom_id < 0)
            error->all(FC_FILE_LINE_FUNC_PARSE, "atom_id have to be non-negative.");

          if ((unsigned) fixed_atom_id >= atom_data->atom_struct_owned.molecule_index.size())
            error->all(FC_FILE_LINE_FUNC_PARSE, "atom_id is larger than the number of atoms");
          std::cout << "fixed_atom_id: " << fixed_atom_id << "\n";
          fixed_molecule_id = atom_data->atom_struct_owned.molecule_index[fixed_atom_id];
          std::cout << "fixed_molecule_id: " << fixed_molecule_id << "\n";

          if (fixed_molecule_id == -1) // adding single atom
          {

            fixed_atoms_id.clear();
            fixed_atoms_id.push_back(fixed_atom_id);

            fixed_atoms_resting_position.clear();
            fixed_atoms_resting_position.push_back(atom_data->atom_struct_owned.position[fixed_atom_id]);
          }
          else // adding the whole molecule
          {
            fixed_atoms_id.clear();
            std::cout << "atom_data->molecule_struct_owned[fixed_molecule_id].atom_list: " << atom_data->molecule_struct_owned[fixed_molecule_id].atom_list.size() << "\n";

            for (auto & i : atom_data->molecule_struct_owned[fixed_molecule_id].atom_list)
            {
              fixed_atoms_id.push_back(i);
            }

            fixed_atoms_resting_position.clear();
            for (auto & i : fixed_atoms_id)
            {
              fixed_atoms_resting_position.push_back(atom_data->atom_struct_owned.position[i]);
            }
          }
          std::cout << "fixed_atoms_id: ";
          for (unsigned int i = 0; i <fixed_atoms_id.size(); ++i)
            std::cout << fixed_atoms_id[i] << " , ";
          std::cout << "\n";
        }
        else if (string_cmp(t, "moving_atom_id"))
        {
          GET_OR_CHOOSE_A_INT(moving_atom_id, "", "")
          if (moving_atom_id < 0)
            error->all(FC_FILE_LINE_FUNC_PARSE, "atom_id have to be non-negative.");

          if ((unsigned) moving_atom_id >= atom_data->atom_struct_owned.molecule_index.size())
            error->all(FC_FILE_LINE_FUNC_PARSE, "atom_id is larger than the number of atoms");

          moving_molecule_id = atom_data->atom_struct_owned.molecule_index[moving_atom_id];
          if (moving_molecule_id == -1) // adding single atom
          {
            moving_atoms_id.clear();
            moving_atoms_id.push_back(moving_atom_id);

            moving_atoms_resting_position.clear();
            moving_atoms_resting_position.push_back(atom_data->atom_struct_owned.position[moving_atom_id]);

            moving_atoms_initial_position.clear();
            moving_atoms_initial_position.push_back(atom_data->atom_struct_owned.position[moving_atom_id]);
          }
          else // adding the whole molecule
          {
            moving_atoms_id.clear();
            for (auto & i : atom_data->molecule_struct_owned[moving_molecule_id].atom_list)
            {
              moving_atoms_id.push_back(i);
            }

            moving_atoms_resting_position.clear();
            moving_atoms_initial_position.clear();

            for (auto & i : moving_atoms_id)
            {
              moving_atoms_resting_position.push_back(atom_data->atom_struct_owned.position[i]);
              moving_atoms_initial_position.push_back(atom_data->atom_struct_owned.position[i]);
            }
          }
          std::cout << "moving_atoms_id: ";
          for (unsigned int i = 0; i <moving_atoms_id.size(); ++i)
            std::cout << moving_atoms_id[i] << " , ";
          std::cout << "\n";
        }
        else if (string_cmp(t, "step"))
        {
          GET_OR_CHOOSE_A_INT(step, "", "")
          if (step <= 0)
            error->all(FC_FILE_LINE_FUNC_PARSE, "step have to be a positive number.");
        }
        else if (string_cmp(t, "set_position"))
        {
          double x = 0, y = 0, z = 0;
          GET_OR_CHOOSE_A_REAL(x, "", "")
          GET_OR_CHOOSE_A_REAL(y, "", "")
          GET_OR_CHOOSE_A_REAL(z, "", "")
          position = Vector3d<double>{x, y, z};
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
        else if (string_cmp(t, "initial_pulling"))
        {
          GET_OR_CHOOSE_A_REAL(initial_pulling, "", "")          
        }
        else if (string_cmp(t, "final_pulling"))
        {
          GET_OR_CHOOSE_A_REAL(final_pulling, "", "")          
        }
        else if (string_cmp(t, "pulling_rate"))
        {
          GET_OR_CHOOSE_A_REAL(pulling_rate, "", "")          
        }
        else if (string_cmp(t, "reaction_coordinate"))
        {
          std::string s = "x";
          GET_A_STRING(s, "", "");
          if (s == "x")
            reaction_coordinate = 'x';
          else if (s == "y")
            reaction_coordinate = 'y';
          else if (s == "z")
            reaction_coordinate = 'z';
          else
            error->all(FC_FILE_LINE_FUNC_PARSE, "Expected x, y or z as reaction coordinate");
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
        else if (string_cmp(t, "open_production_files"))
        {
          open_production_files();
        }
        else if (string_cmp(t, "close_production_files"))
        {
          close_production_files();
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
        // else if (string_cmp(t, "set_writer_atom_data"))
        // {
        //   FIND_OBJECT_BY_NAME(writer, it)
        //   FC_CHECK_OBJECT_CLASS_NAME(writer, it, atom_data)
        //   auto a = *dynamic_cast<writer::Atom_data *>(object_container->writer[it->second.index]);         
          
        // }
        else
          FC_ERR_UNDEFINED_VAR(t)
      }

      return in_file;
    }

    void Umbrella_sampling_g::fix_file_prefix()
    {
      if (!is_absolute_path(file_prefix))
      {
        if (file_prefix == "")
          file_prefix = get_current_directory();
        else
          file_prefix = join_path(get_current_directory(), file_prefix);
      }
    }


    void Umbrella_sampling_g::init_metadata()
    {
      metadata_mode = true;


      if (writerXYZ == nullptr)
        writerXYZ = new writer::Atom_data(caviar_);

      writerXYZ->output_xyz = true;
      writerXYZ->atom_data = atom_data;
      writerXYZ->domain = domain;
      writerXYZ->xyz_step = 1;
      mkdir(file_prefix.c_str(), 0777);
      ofs_summary_distance.open(join_path(file_prefix , file_name_summary_distance));
      ofs_summary_distance << "Frame" <<'\t' << "Distance" << '\n';

      pulling_value = initial_pulling;
    }

    void Umbrella_sampling_g::finish_metadata()
    {
      metadata_mode = false;
      if (writerXYZ != nullptr)
        delete writerXYZ;
      writerXYZ = nullptr;


      ofs_summary_distance.close();
    }

    void Umbrella_sampling_g::metadata_function()
    {
      if (!metadata_mode)
        return;      



      //int i = atom_data->atom_id_to_index[atom_id];
      //auto p_i = atom_data->atom_struct_owned.position[i];

      ofs_data << step_counter;
      
      if (pulling_value < final_pulling)
      {
        pulling_value += pulling_rate;
      }
      else
      {
        pulling_value = final_pulling;
      }

      for (unsigned int f_i = 0; f_i < moving_atoms_id.size(); ++f_i)
      {
        //int i = atom_data->atom_id_to_index[moving_atoms_id[f_i]];        

        if (reaction_coordinate == 'x')
        {
          moving_atoms_resting_position[f_i].x = moving_atoms_initial_position[f_i].x + pulling_value;
        }
        else if (reaction_coordinate == 'y')
        {
          moving_atoms_resting_position[f_i].y = moving_atoms_initial_position[f_i].y + pulling_value;
        }
        else if (reaction_coordinate == 'z')
        {
          moving_atoms_resting_position[f_i].z = moving_atoms_initial_position[f_i].z + pulling_value;
        }
        else
        {
          error->all(FC_FILE_LINE_FUNC, "Expected x, y or z as reaction coordinate");
        }
      }

      if (step_counter % step == 0)
      {
        writerXYZ->file_name_xyz =  join_path(file_prefix , "conf_" + std::to_string(step_counter));
        writerXYZ->initialize();
        writerXYZ->write(step_counter,1);
        writerXYZ->close_files();

        auto com_fixed = atom_data->owned_position_cm(fixed_atoms_id);
        auto com_moving = atom_data->owned_position_cm(moving_atoms_id);

        auto com_diff = com_fixed - com_moving;

        ofs_summary_distance << step_counter <<'\t' << norm(com_diff) << '\n';
      }
        
    }
    
    void Umbrella_sampling_g::open_production_files()
    {
      
      mkdir(file_prefix.c_str(), 0777);
      std::string file_name_tmp = join_path(file_prefix, file_name_metadata + ".txt");
      ofs_metadata.open(file_name_tmp.c_str());

      std::string file_name_stat_tmp = join_path(file_prefix, file_name_stat + ".txt");
      ofs_stat.open(file_name_stat_tmp.c_str());
    }

    void Umbrella_sampling_g::close_production_files()
    {
      ofs_metadata.close();
      ofs_stat.close();
    }

    void Umbrella_sampling_g::init_production()
    {
		for (size_t i = 0; i < fixed_atoms_id.size(); ++i)
		{
			int id = fixed_atoms_id[i];
			fixed_atoms_resting_position[i] = atom_data->atom_struct_owned.position[id];
		}
		
		for (size_t i = 0; i < moving_atoms_id.size(); ++i)
		{
			int id = moving_atoms_id[i];
			moving_atoms_resting_position[i] = atom_data->atom_struct_owned.position[id];
			moving_atoms_initial_position[i] = atom_data->atom_struct_owned.position[id];
		}
		
      production_mode = true;

      // mkdir(file_prefix.c_str(), 0777);
      // std::string file_name_tmp = join_path(file_prefix, file_name_metadata + ".txt");
      // ofs_metadata.open(file_name_tmp.c_str(), std::ios_base::app);

      // std::string file_name_stat_tmp = join_path(file_prefix, file_name_stat + ".txt");
      // ofs_stat.open(file_name_stat_tmp.c_str(), std::ios_base::app);


      step_counter = -1;

      std::string file_name_tmp = join_path(file_prefix, file_name_data + std::to_string(file_counter) + ".txt");
      ofs_data.open(file_name_tmp.c_str());
      // if (metadata_mode == false)
      // {
      //   error->all(FC_FILE_LINE_FUNC, "Expected metadata activation, i.e., 'init_metadata' call");
      // }


      //------------------------
        //auto com_fixed = atom_data->owned_position_cm(fixed_atoms_id);
        auto com_moving = atom_data->owned_position_cm(moving_atoms_id);
        position = com_moving;
        //auto com_diff = com_fixed - com_moving;

      //------------------------


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
        error->all(FC_FILE_LINE_FUNC, "Expected x, y or z as reaction coordinate");
      }
      ofs_metadata << " " << elastic_coef;
      ofs_metadata << " " << "0"; // coecorrelation time for your time series used in monte-carlo

      if (temperature_is_set)
      {
        ofs_metadata << " " << temperature;
      }
      ofs_metadata << "\n";
    }

    void Umbrella_sampling_g::production_function()
    {
      if (!production_mode)
        return;

      

      if (step_counter % step != 0)
        return;

      // int i = atom_data->atom_id_to_index[atom_id];
      // auto p_i = atom_data->atom_struct_owned.position[i];
      auto p_i = atom_data->owned_position_cm(moving_atoms_id);

      ofs_data << step_counter;

      if (reaction_coordinate == 'x')
      {
        ofs_data << " " << p_i.x;
        if (position_min > p_i.x)
          position_min = p_i.x;
        if (position_max < p_i.x)
          position_max = p_i.x;
        position_mean += p_i.x;
      }
      else if (reaction_coordinate == 'y')
      {
        ofs_data << " " << p_i.y;
        if (position_min > p_i.y)
          position_min = p_i.y;
        if (position_max < p_i.y)
          position_max = p_i.y;
        position_mean += p_i.y;
      }
      else if (reaction_coordinate == 'z')
      {
        ofs_data << " " << p_i.z;
        if (position_min > p_i.z)
          position_min = p_i.z;
        if (position_max < p_i.z)
          position_max = p_i.z;
        position_mean += p_i.z;
      }
      else
      {
        error->all(FC_FILE_LINE_FUNC, "Expected x, y or z as reaction coordinate");
      }
      position_mean_count++;
      ofs_data << "\n";

      //
      // ofs_data << atom_data->atom_struct_owned.position[i] << " " << dr << " " << 0.5 * elastic_coef * dr * dr << "\n";
    }

    void Umbrella_sampling_g::finish_production()
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
      ofs_stat << " " << position_mean << " " << position_min << " " << position_max << "\n"
               << std::flush;

      ofs_data << std::flush;
      production_mode = false;
      ofs_data.close();
      file_counter++;


    }

    void Umbrella_sampling_g::verify_settings()
    {
      FC_NULLPTR_CHECK(atom_data)
      FC_NULLPTR_CHECK(domain)
      // my_mpi_rank = atom_data->get_mpi_rank();
      // if ((unsigned int)atom_id > atom_data->atom_id_to_index.size() - 1)
      //   error->all(FC_FILE_LINE_FUNC, "Invalid atom_id:" + std::to_string(atom_id));

      // int i = atom_data->atom_id_to_index[atom_id];
    }

    void Umbrella_sampling_g::calculate_acceleration()
    {
      FC_OBJECT_VERIFY_SETTINGS

      step_counter++;

      auto &type = atom_data->atom_struct_owned.type;
      auto &mass_inv = atom_data->atom_type_params.mass_inv;

      for (unsigned int f_i = 0; f_i < fixed_atoms_id.size(); ++f_i)
      {
        int i = atom_data->atom_id_to_index[fixed_atoms_id[f_i]];

#if defined(CAVIAR_WITH_MPI)
        dr = fixed_atoms_resting_position[f_i] - atom_data->atom_struct_owned.position[i];
#else
        dr = domain->periodic_distance(fixed_atoms_resting_position[f_i] - atom_data->atom_struct_owned.position[i]);
#endif
        const auto force = -elastic_coef * dr;

        atom_data->atom_struct_owned.acceleration[i] -= force * mass_inv[type[i]];
      }

      for (unsigned int f_i = 0; f_i < moving_atoms_id.size(); ++f_i)
      {
        int i = atom_data->atom_id_to_index[moving_atoms_id[f_i]];

#if defined(CAVIAR_WITH_MPI)
        dr = moving_atoms_resting_position[f_i] - atom_data->atom_struct_owned.position[i];
#else
        dr = domain->periodic_distance(moving_atoms_resting_position[f_i] - atom_data->atom_struct_owned.position[i]);
#endif
        const auto force = -elastic_coef * dr;

        atom_data->atom_struct_owned.acceleration[i] -= force * mass_inv[type[i]];
      }

      metadata_function();
      production_function();
    }

  } // force_field

}
