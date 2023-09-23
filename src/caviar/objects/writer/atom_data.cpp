
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
#include "caviar/interpreter/communicator.h"

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

    if (ofs_xyz_mpi.is_open())
      ofs_xyz_mpi.close();

    if (ofs_xyz_ghost.is_open())
      ofs_xyz_ghost.close();

    if (ofs_xyz_ghost_mpi.is_open())
      ofs_xyz_ghost_mpi.close();

    if (ofs_energy.is_open())
      ofs_energy.close();

    if (ofs_energy_mpi.is_open())
      ofs_energy_mpi.close();

    if (ofs_temperature.is_open())
      ofs_temperature.close();

    if (ofs_temperature_mpi.is_open())
      ofs_temperature_mpi.close();

    if (ofs_pressure.is_open())
      ofs_pressure.close();

    if (ofs_pressure_mpi.is_open())
      ofs_pressure_mpi.close();

    if (ofs_povray.is_open())
      ofs_povray.close();

    if (ofs_povray_mpi.is_open())
      ofs_povray_mpi.close();

    if (ofs_msd.is_open())
      ofs_msd.close();

    if (ofs_msd_mpi.is_open())
      ofs_msd_mpi.close();

    if (ofs_volume.is_open())
      ofs_volume.close();

    if (ofs_volume_mpi.is_open())
      ofs_volume_mpi.close();
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
      else if (string_cmp(t, "xyz_mpi_per_process"))
      {
        GET_OR_CHOOSE_A_INT(xyz_mpi_per_process, "", "")
      }
      else if (string_cmp(t, "xyz_ghost_mpi_per_process"))
      {
        GET_OR_CHOOSE_A_INT(xyz_ghost_mpi_per_process, "", "")
      }
      else if (string_cmp(t, "energy_mpi_per_process"))
      {
        GET_OR_CHOOSE_A_INT(energy_mpi_per_process, "", "")
      }
      else if (string_cmp(t, "temperature_mpi_per_process"))
      {
        GET_OR_CHOOSE_A_INT(temperature_mpi_per_process, "", "")
      }
      else if (string_cmp(t, "pressure_mpi_per_process"))
      {
        GET_OR_CHOOSE_A_INT(pressure_mpi_per_process, "", "")
      }
      else if (string_cmp(t, "povray_mpi_per_process"))
      {
        GET_OR_CHOOSE_A_INT(povray_mpi_per_process, "", "")
      }
      else if (string_cmp(t, "msd_mpi_per_process"))
      {
        GET_OR_CHOOSE_A_INT(msd_mpi_per_process, "", "")
      }
      else if (string_cmp(t, "volume_mpi_per_process"))
      {
        GET_OR_CHOOSE_A_INT(volume_mpi_per_process, "", "")
      }
      else if (string_cmp(t, "xyz_mpi_rank0"))
      {
        GET_OR_CHOOSE_A_INT(xyz_mpi_rank0, "", "")
      }
      else if (string_cmp(t, "xyz_ghost_mpi_rank0"))
      {
        GET_OR_CHOOSE_A_INT(xyz_ghost_mpi_rank0, "", "")
      }
      else if (string_cmp(t, "energy_mpi_rank0"))
      {
        GET_OR_CHOOSE_A_INT(energy_mpi_rank0, "", "")
      }
      else if (string_cmp(t, "temperature_mpi_rank0"))
      {
        GET_OR_CHOOSE_A_INT(temperature_mpi_rank0, "", "")
      }
      else if (string_cmp(t, "pressure_mpi_rank0"))
      {
        GET_OR_CHOOSE_A_INT(pressure_mpi_rank0, "", "")
      }
      else if (string_cmp(t, "povray_mpi_rank0"))
      {
        GET_OR_CHOOSE_A_INT(povray_mpi_rank0, "", "")
      }
      else if (string_cmp(t, "msd_mpi_rank0"))
      {
        GET_OR_CHOOSE_A_INT(msd_mpi_rank0, "", "")
      }
      else if (string_cmp(t, "volume_mpi_rank0"))
      {
        GET_OR_CHOOSE_A_INT(volume_mpi_rank0, "", "")
      }
      else if (string_cmp(t, "xyz_ghost_step"))
      {
        GET_OR_CHOOSE_A_INT(xyz_ghost_step, "", "")
        output_xyz_ghost = true;
        if (xyz_ghost_step <= 0)
        {
          output_xyz_ghost = false;
        }
      }
      else if (string_cmp(t, "xyz_step"))
      {
        GET_OR_CHOOSE_A_INT(xyz_step, "", "")
        output_xyz = true;
        if (xyz_step <= 0)
        {
          output_xyz = false;
        }
      }
      else if (string_cmp(t, "povray_step"))
      {
        GET_OR_CHOOSE_A_INT(povray_step, "", "")
        output_povray = true;
        if (povray_step <= 0)
        {
          output_povray = false;
        }
      }
      else if (string_cmp(t, "energy_step"))
      {
        GET_OR_CHOOSE_A_INT(energy_step, "", "")
        output_energy = true;
        if (energy_step <= 0)
        {
          output_energy = false;
        }
      }
      else if (string_cmp(t, "temperature_step"))
      {
        GET_OR_CHOOSE_A_INT(temperature_step, "", "")
        output_temperature = true;
        if (temperature_step <= 0)
        {
          output_temperature = false;
        }
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
        if (msd_step <= 0)
        {
          output_msd = false;
        }
      }
      else if (string_cmp(t, "msd_initial_step"))
      {
        GET_OR_CHOOSE_A_INT(msd_initial_step, "", "")
        // output_msd = true;
      }
      else if (string_cmp(t, "volume_step"))
      {
        GET_OR_CHOOSE_A_INT(volume_step, "", "")
        output_volume = true;
        if (pressure_step <= 0)
        {
          output_pressure = false;
        }
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
      else if (string_cmp(t, "file_name_volume"))
      {
        file_name_volume = parser->get_val_token().string_value;
      }
      else if (string_cmp(t, "open_files"))
      {
        open_files();
      }
      else if (string_cmp(t, "close_files"))
      {
        close_files();
      }
      else if (string_cmp(t, "xyz_output_id"))
      {
        GET_OR_CHOOSE_A_INT(xyz_output_id, "", "")
      }
      else if (string_cmp(t, "xyz_output_velocity"))
      {
        GET_OR_CHOOSE_A_INT(xyz_output_velocity, "", "")
      }
      else if (string_cmp(t, "xyz_output_acceleration"))
      {
        GET_OR_CHOOSE_A_INT(xyz_output_acceleration, "", "")

        // std::ofstream ofs("o_acc");
        // const auto &pos = atom_data->atom_struct_owned.position;
        // const auto &acc = atom_data->atom_struct_owned.acceleration;
        // for (unsigned int i = 0; i < pos.size(); ++i)
        // {
        //   ofs << i << " " << acc[i].x << "\t" << acc[i].y << "\t" << acc[i].z << "\n";
        // }
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

    if (output_volume)
    {
      FC_NULLPTR_CHECK(domain)
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

      if (output_volume)
        if (!ofs_volume.is_open())
          ofs_volume.open((file_name_volume + ".txt").c_str());
    }
#elif defined(CAVIAR_WITH_MPI)

    if (output_xyz && xyz_mpi_per_process)
      if (!ofs_xyz_mpi.is_open())
        ofs_xyz_mpi.open((file_name_xyz + "_mpi" + std::to_string(my_mpi_rank) + ".xyz").c_str());

    if (output_xyz_ghost && xyz_ghost_mpi_per_process)
      if (!ofs_xyz_ghost_mpi.is_open())
        ofs_xyz_ghost_mpi.open((file_name_xyz_ghost + "_mpi" + std::to_string(my_mpi_rank) + ".xyz").c_str());

    if (output_povray && povray_mpi_per_process)
      if (!ofs_povray_mpi.is_open())
        ofs_povray_mpi.open((file_name_povray + "_mpi" + std::to_string(my_mpi_rank) + ".pov").c_str());

    if (output_energy && energy_mpi_per_process)
      if (!ofs_energy_mpi.is_open())
        ofs_energy_mpi.open((file_name_energy + "_mpi" + std::to_string(my_mpi_rank) + ".txt").c_str());

    if (output_temperature && temperature_mpi_per_process)
      if (!ofs_temperature_mpi.is_open())
        ofs_temperature_mpi.open((file_name_temperature + "_mpi" + std::to_string(my_mpi_rank) + ".txt").c_str());

    if (output_pressure && pressure_mpi_per_process)
      if (!ofs_pressure_mpi.is_open())
        ofs_pressure_mpi.open((file_name_pressure + "_mpi" + std::to_string(my_mpi_rank) + ".txt").c_str());

    if (output_msd && msd_mpi_per_process)
      if (!ofs_msd_mpi.is_open())
        ofs_msd_mpi.open((file_name_msd + "_mpi" + std::to_string(my_mpi_rank) + ".txt").c_str());

    if (output_volume && volume_mpi_per_process)
      if (!ofs_volume_mpi.is_open())
        ofs_volume_mpi.open((file_name_volume + "_mpi" + std::to_string(my_mpi_rank) + ".txt").c_str());

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

      if (output_volume)
        if (!ofs_volume.is_open())
          ofs_volume.open((file_name_volume + ".txt").c_str());
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

    if (output_volume)
      if (!ofs_volume.is_open())
        ofs_volume.open((file_name_volume + ".txt").c_str());
#endif
  }

  void Atom_data::open_files() {}
  void Atom_data::close_files() {}
  void Atom_data::generate() {}

  void Atom_data::report_xyz_dump(int64_t i, double)
  {
    if (my_mpi_rank != 0)
      return;

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
      write_mpi_shared_atoms(i, t);
    }

    write_serial(i, t);

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

    if (output_volume)
      if (i % volume_step == 0)
        dump_volume_mpi(i, t);
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

    if (output_volume)
      if (i % volume_step == 0)
        dump_volume_mpi_shared_atoms(i, t);
  }

  void Atom_data::write_serial(int64_t i, double t)
  {
    if (output_xyz)
      if (i % xyz_step == 0)
      {
        dump_xyz_serial(i, t);
        report_xyz_dump(i, t);
      }

    if (output_xyz_ghost)
      if (i % xyz_ghost_step == 0)
        dump_xyz_ghost_serial(i, t);

    if (output_energy)
      if (i % energy_step == 0)
        dump_energy_serial(i, t);

    if (output_temperature)
      if (i % temperature_step == 0)
        dump_temperature_serial(i, t);

    if (output_pressure)
      if (i % pressure_step == 0)
        dump_pressure_serial(i, t);

    if (output_povray)
      if (i % povray_step == 0)
        dump_povray_serial(i, t);

    if (output_msd)
      if (i % msd_step == 0)
        dump_msd_serial(i, t);

    if (output_volume)
      if (i % volume_step == 0)
        dump_volume_serial(i, t);
  }

  void Atom_data::start_new_files() {}              // add_time_to_previous
  void Atom_data::start_new_files(std::string &) {} // add_time_to_previous

} // Atom_data

CAVIAR_NAMESPACE_CLOSE
