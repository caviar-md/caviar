
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

#include "caviar/objects/writer/atom_data.h"
#include "caviar/objects/atom_data.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/utility/time_utility.h"
#include "caviar/objects/unique/time_function_3d.h"
#include "caviar/objects/domain.h"

#include <ctime>
#include <sys/stat.h> // used for mkdir()

CAVIAR_NAMESPACE_OPEN

namespace writer
{

  Atom_data::Atom_data(CAVIAR *fptr) : Writer{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    wallTimeXyzDump1 = get_wall_time();
  }

  Atom_data::~Atom_data()
  {
    if (ofs_xyz.is_open())
      ofs_xyz.close();

    if (ofs_xyz_ghost.is_open())
      ofs_xyz_ghost.close();

    if (ofs_xyz_mpi.is_open())
      ofs_xyz_mpi.close();

    if (ofs_xyz_ghost_mpi.is_open())
      ofs_xyz_ghost_mpi.close();

    if (ofs_energy.is_open())
      ofs_energy.close();

    if (ofs_temperature.is_open())
      ofs_temperature.close();

    if (ofs_pressure.is_open())
      ofs_pressure.close();

    if (ofs_povray.is_open())
      ofs_povray.close();

    if (ofs_msd.is_open())
      ofs_msd.close();
  }

  bool Atom_data::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else if (string_cmp(t, "set_domain") || string_cmp(t, "domain"))
      {
        FIND_OBJECT_BY_NAME(domain, it)
        domain = object_container->domain[it->second.index];
      }
      else if (string_cmp(t, "mpi_separate_files"))
      {
        mpi_separate_files = true;
      }
      else if (string_cmp(t, "mpi_single_file"))
      {
        mpi_single_file = true;
      }
      else if (string_cmp(t, "xyz_ghost_step"))
      {
        GET_OR_CHOOSE_A_INT(xyz_ghost_step, "", "")
        output_xyz_ghost = true;
      }
      else if (string_cmp(t, "xyz_step"))
      {
        GET_OR_CHOOSE_A_INT(xyz_step, "", "")
        output_xyz = true;
        xyz_ghost_step = xyz_step;
      }
      else if (string_cmp(t, "povray_step"))
      {
        GET_OR_CHOOSE_A_INT(povray_step, "", "")
        output_povray = true;
      }
      else if (string_cmp(t, "energy_step"))
      {
        GET_OR_CHOOSE_A_INT(energy_step, "", "")
        output_energy = true;
      }
      else if (string_cmp(t, "temperature_step"))
      {
        GET_OR_CHOOSE_A_INT(temperature_step, "", "")
        output_temperature = true;
      }
      else if (string_cmp(t, "pressure_step"))
      {
        GET_OR_CHOOSE_A_INT(pressure_step, "", "")
        output_pressure = true;
        if (pressure_step <= 0)
        {
          output_pressure = false;
        }
      }
      else if (string_cmp(t, "msd_step"))
      {
        GET_OR_CHOOSE_A_INT(msd_step, "", "")
        output_msd = true;
      }
      else if (string_cmp(t, "msd_initial_step"))
      {
        GET_OR_CHOOSE_A_INT(msd_initial_step, "", "")
        output_msd = true;
      }
      else if (string_cmp(t, "file_name_xyz"))
      {
        file_name_xyz = parser->get_val_token().string_value;
      }
      else if (string_cmp(t, "file_name_xyz_ghost"))
      {
        file_name_xyz_ghost = parser->get_val_token().string_value;
      }
      else if (string_cmp(t, "file_name_energy"))
      {
        file_name_energy = parser->get_val_token().string_value;
      }
      else if (string_cmp(t, "file_name_temperature"))
      {
        file_name_temperature = parser->get_val_token().string_value;
      }
      else if (string_cmp(t, "file_name_pressure"))
      {
        file_name_pressure = parser->get_val_token().string_value;
      }
      else if (string_cmp(t, "file_name_povray"))
      {
        file_name_povray = parser->get_val_token().string_value;
      }
      else if (string_cmp(t, "file_name_msd"))
      {
        file_name_msd = parser->get_val_token().string_value;
      }
      else if (string_cmp(t, "open_files"))
      {
        open_files();
      }
      else if (string_cmp(t, "close_files"))
      {
        close_files();
      }
      else if (string_cmp(t, "output_id"))
      {
        output_id = true;
      }
      else if (string_cmp(t, "output_velocity"))
      {
        output_velocity = true;
      }
      else if (string_cmp(t, "output_acceleration"))
      {
        // output_acceleration = true;
        std::ofstream ofs("o_acc");
        const auto &pos = atom_data->atom_struct_owned.position;
        const auto &acc = atom_data->atom_struct_owned.acceleration;
        for (unsigned int i = 0; i < pos.size(); ++i)
        {
          ofs << i << " " << acc[i].x << "\t" << acc[i].y << "\t" << acc[i].z << "\n";
        }
      }
      else if (string_cmp(t, "set_position_offset"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        FC_CHECK_OBJECT_CLASS_NAME(unique, it, time_function_3d)
        unique::Time_function_3d *a = dynamic_cast<unique::Time_function_3d *>(object_container->unique[it->second.index]);
        position_offset = a;
      }
      else
        FC_ERR_UNDEFINED_VAR(t)
    }
    return in_file;
  }

  void Atom_data::initialize()
  {
    initialized = true;

    FC_NULLPTR_CHECK(atom_data)

    if (output_msd)
    {
      FC_NULLPTR_CHECK(domain)
      if (!atom_data->get_msd_process())
        error->all(FC_FILE_LINE_FUNC, "In order to have 'output_msd' in writer::Atom_data, 'msd_process' must be activated in atom_data::Atom_data");
    }

    // --- just to make povray outpuy folder ---
    if (my_mpi_rank == 0)
    {
      if (output_povray)
      {
        std::string str_folder_pov;
        str_folder_pov.append("o_pov");
        const char *char_folder_pov = str_folder_pov.c_str();
        mkdir(char_folder_pov, 0777); // make povray  output folder //
      }
    }
    /*
    #ifdef CAVIAR_WITH_MPI
      sprintf ( buffer, "_me%u", comm->me );
      str_filename.append ( buffer);
    #endif
    */

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
    if (my_mpi_rank == 0)
    {
      if (output_xyz)
        if (!ofs_xyz.is_open())
          ofs_xyz.open((file_name_xyz + ".xyz").c_str());

      if (output_xyz_ghost)
        if (!ofs_xyz_ghost.is_open())
          ofs_xyz_ghost.open((file_name_xyz_ghost + ".xyz").c_str());

      if (output_povray)
        if (!ofs_povray.is_open())
          ofs_povray.open((file_name_povray + ".pov").c_str());

      if (output_energy)
        if (!ofs_energy.is_open())
          ofs_energy.open((file_name_energy + ".txt").c_str());

      if (output_temperature)
        if (!ofs_temperature.is_open())
          ofs_temperature.open((file_name_temperature + ".txt").c_str());

      if (output_pressure)
        if (!ofs_pressure.is_open())
          ofs_pressure.open((file_name_pressure + ".txt").c_str());

      if (output_msd)
        if (!ofs_msd.is_open())
          ofs_msd.open((file_name_msd + ".txt").c_str());
    }
#elif defined(CAVIAR_WITH_MPI)
    if (mpi_separate_files)
    {
      if (output_xyz)
        if (!ofs_xyz_mpi.is_open())
          ofs_xyz_mpi.open((file_name_xyz + "_mpi" + std::to_string(my_mpi_rank) + ".xyz").c_str());

      if (output_xyz_ghost)
        if (!ofs_xyz_ghost_mpi.is_open())
          ofs_xyz_ghost_mpi.open((file_name_xyz_ghost + "_mpi" + std::to_string(my_mpi_rank) + ".xyz").c_str());
    }
    if (my_mpi_rank == 0)
    {
      if (mpi_single_file)
      {
        if (output_xyz)
          if (!ofs_xyz.is_open())
            ofs_xyz.open((file_name_xyz + ".xyz").c_str());

        if (output_xyz_ghost)
          if (!ofs_xyz_ghost.is_open())
            ofs_xyz_ghost.open((file_name_xyz_ghost + ".xyz").c_str());
      }
      if (output_povray)
        if (!ofs_povray.is_open())
          ofs_povray.open((file_name_povray + ".pov").c_str());

      if (output_energy)
        if (!ofs_energy.is_open())
          ofs_energy.open((file_name_energy + ".txt").c_str());

      if (output_temperature)
        if (!ofs_temperature.is_open())
          ofs_temperature.open((file_name_temperature + ".txt").c_str());

      if (output_pressure)
        if (!ofs_pressure.is_open())
          ofs_pressure.open((file_name_pressure + ".txt").c_str());

      if (output_msd)
        if (!ofs_msd.is_open())
          ofs_msd.open((file_name_msd + ".txt").c_str());
    } 
#else

    if (output_xyz)
      if (!ofs_xyz.is_open())
        ofs_xyz.open((file_name_xyz + ".xyz").c_str());

    if (output_xyz_ghost)
      if (!ofs_xyz_ghost.is_open())
        ofs_xyz_ghost.open((file_name_xyz_ghost + ".xyz").c_str());

    if (output_povray)
      if (!ofs_povray.is_open())
        ofs_povray.open((file_name_povray + ".pov").c_str());

    if (output_energy)
      if (!ofs_energy.is_open())
        ofs_energy.open((file_name_energy + ".txt").c_str());

    if (output_temperature)
      if (!ofs_temperature.is_open())
        ofs_temperature.open((file_name_temperature + ".txt").c_str());

    if (output_pressure)
      if (!ofs_pressure.is_open())
        ofs_pressure.open((file_name_pressure + ".txt").c_str());

    if (output_msd)
      if (!ofs_msd.is_open())
        ofs_msd.open((file_name_msd + ".txt").c_str());
#endif
  }

  void Atom_data::open_files() {}
  void Atom_data::close_files() {}
  void Atom_data::generate() {}

  void Atom_data::report_xyz_dump(int64_t i, double t)
  {
    double wallTimeXyzDump2 = get_wall_time();

    double dtstart = wallTimeXyzDump2 - wallTimeXyzDump1;
    std::string s = "dump_xyz [" + std::to_string(i) +
                    +"] (" + std::to_string(dtstart) + " S)";
    output->info(s, 2);
    wallTimeXyzDump1 = wallTimeXyzDump2;
  }

  void Atom_data::write(int64_t i, double t)
  {

    if (!initialized)
      initialize();


#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

    if (my_mpi_rank == 0)
    {
      write_serial(i, t);
    }
#elif defined(CAVIAR_WITH_MPI)

    if (comm->nprocs > 1)
    {
      if (mpi_single_file)
      {
        write_mpi_shared_atoms(i, t);
        // write_mpi(i);
      }
      if (mpi_separate_files)
        write_serial(i, t, true);
    }
    else
    {
      write_serial(i, t);
    }

#else

    write_serial(i, t);

#endif




  }

  void Atom_data::write_mpi(int64_t i, double t)
  {
    if (output_xyz)
      if (i % xyz_step == 0)
      {
        dump_xyz_mpi(i, t);
        report_xyz_dump(i, t);
      }

    if (output_xyz_ghost)
      if (i % xyz_ghost_step == 0)
        dump_xyz_ghost_mpi(i, t);

    if (output_energy)
      if (i % energy_step == 0)
        dump_energy_mpi(i, t);

    if (output_temperature)
      if (i % temperature_step == 0)
        dump_temperature_mpi(i, t);

    if (output_pressure)
      if (i % pressure_step == 0)
        dump_pressure_mpi(i, t);

    if (output_povray)
      if (i % povray_step == 0)
        dump_povray_mpi(i, t);

    if (output_msd)
      if (i % msd_step == 0)
        dump_msd_mpi(i, t);
  }

  void Atom_data::write_mpi_shared_atoms(int64_t i, double t)
  {
    if (output_xyz)
      if (i % xyz_step == 0)
      {
        dump_xyz_mpi_shared_atoms(i, t);
        report_xyz_dump(i, t);
      }

    if (output_xyz_ghost)
      if (i % xyz_ghost_step == 0)
        dump_xyz_ghost_mpi_shared_atoms(i, t);

    if (output_energy)
      if (i % energy_step == 0)
        dump_energy_mpi_shared_atoms(i, t);

    if (output_temperature)
      if (i % temperature_step == 0)
        dump_temperature_mpi_shared_atoms(i, t);

    if (output_pressure)
      if (i % pressure_step == 0)
        dump_pressure_mpi_shared_atoms(i, t);

    if (output_povray)
      if (i % povray_step == 0)
        dump_povray_mpi_shared_atoms(i, t);

    if (output_msd)
      if (i % msd_step == 0)
        dump_msd_mpi_shared_atoms(i, t);
  }

  void Atom_data::write_serial(int64_t i, double t, bool mpi_files)
  {
    if (output_xyz)
      if (i % xyz_step == 0)
      {
        dump_xyz_serial(i, t, mpi_files);
        report_xyz_dump(i, t);
      }

    if (output_xyz_ghost)
      if (i % xyz_ghost_step == 0)
        dump_xyz_ghost_serial(i, t, mpi_files);

    if (output_energy)
      if (i % energy_step == 0)
        dump_energy_serial(i, t, mpi_files);

    if (output_temperature)
      if (i % temperature_step == 0)
        dump_temperature_serial(i, t, mpi_files);

    if (output_pressure)
      if (i % pressure_step == 0)
        dump_pressure_serial(i, t, mpi_files);

    if (output_povray)
      if (i % povray_step == 0)
        dump_povray_serial(i, t, mpi_files);

    if (output_msd)
      if (i % msd_step == 0)
        dump_msd_serial(i, t, mpi_files);
  }

  void Atom_data::start_new_files() {}              // add_time_to_previous
  void Atom_data::start_new_files(std::string &) {} // add_time_to_previous

} // Atom_data

CAVIAR_NAMESPACE_CLOSE
