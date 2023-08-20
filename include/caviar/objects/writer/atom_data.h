
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

#ifndef CAVIAR_OBJECTS_WRITER_ATOMDATA_H
#define CAVIAR_OBJECTS_WRITER_ATOMDATA_H

#include "caviar/objects/writer.h"
#include "caviar/objects/atom_data/utility/mpi_packet_info.h"

CAVIAR_NAMESPACE_OPEN

class Atom_data;
class Domain;
namespace unique
{
  class Time_function_3d;
}
namespace writer
{

  /**
   * This class has a writer for Atom_data class
   *
   *
   */
  class Atom_data : public Writer
  {
  public:
    Atom_data(class CAVIAR *);
    ~Atom_data();
    bool read(class caviar::interpreter::Parser *);
    void initialize();
    void write();
    void write(int64_t);                 // current time_step
    void write(double);                  // current time
    void write(int64_t, double);         // time_step and time
    void start_new_files();              // add_time_to_previous
    void start_new_files(std::string &); // add_time_to_previous
    void open_files();
    void close_files();
    void generate();
    void dump_kinetic_energy();

    unique::Time_function_3d *position_offset = nullptr;

    /////////////

    caviar::Atom_data *atom_data = nullptr;

    // used in msd calculations
    Domain *domain = nullptr;

    void dump_energy(int64_t);         // dump energy to file
    void dump_energy(int64_t, double); // dump energy to file

    void dump_temperature(int64_t);         // dump temperature to file
    void dump_temperature(int64_t, double); // dump temperature to file


    void dump_xyz(int64_t);         // dump positions to file in xyz format
    void dump_xyz(int64_t, double); // dump positions to file in xyz format

    void dump_xyz_serial(int64_t, bool mpi_files = false);   // dump positions to file in xyz format.
    void dump_xyz_mpi(int64_t);                            // dump positions to file in xyz format. Number of atoms are different in mpi domains
    void dump_xyz_mpi_shared_atoms(int64_t);               // dump positions to file in xyz format. Number of atoms are the same in mpi domains

    void dump_xyz_ghost(int64_t);         // dump positions to file in xyz format

    void dump_xyz_ghost_serial(int64_t, bool mpi_files = false);      // dump positions to file in xyz format
    void dump_xyz_ghost_mpi(int64_t);                               // dump positions to file in xyz format. Number of atoms are different in mpi domains

    void dump_povray(int64_t);         // dump positions to snapshot files in povray format
    void dump_povray(int64_t, double); // dump positions to snapshot files in povray format

    void dump_msd(int64_t);         //
    void dump_msd(int64_t, double); //

    int64_t energy_step = 100;
    int64_t temperature_step = 100;
    int64_t xyz_step = 100;
    int64_t xyz_ghost_step = 100;
    int64_t povray_step = 100;
    int64_t msd_step = 100; // number of steps to output data

    int msd_type; // type of atom that should be used in msd calculations.
    int msd_initial_step = 0;

    std::ofstream ofs_energy;
    std::ofstream ofs_temperature;
    std::ofstream ofs_xyz;
    std::ofstream ofs_xyz_ghost;
    std::ofstream ofs_xyz_mpi;
    std::ofstream ofs_xyz_ghost_mpi;
    std::ofstream ofs_velocities;
    std::ofstream ofs_povray;
    std::ofstream ofs_msd; // mean square distance

    // if true, outputs would be created    
    bool output_energy = false;
    bool output_temperature = false;
    bool output_xyz = false;
    bool output_xyz_ghost = false;
    bool output_povray = false;
    bool output_msd = false;

    std::string file_name_xyz = "o_xyz";
    std::string file_name_xyz_ghost = "o_xyz_g";
    std::string file_name_energy = "o_energy";
    std::string file_name_temperature = "o_temperature";
    std::string file_name_povray = "o_pov";
    std::string file_name_msd = "o_msd";

    // dump velocity and acceleration alongside position in xyz file
    bool output_velocity = false; //xyz
    bool output_acceleration = false; //xyz
    bool output_id = false; //xyz

    bool mpi_separate_files = false;
    bool mpi_single_file = false;

    std::vector<caviar::Vector<double>> msd_initial_position;

    // records previous wallTime of XYZ dump.
    double wallTimeXyzDump1;

    /////////

  /**
   * It stores the info about location of the owned Atom_struct data inside the MPI packets.
   */
  atom_data::MPI_packet_info mpi_packet_info_owned;

  /**
   * It stores the info about location of the ghost Atom_struct data inside the MPI packets.
   */
  atom_data::MPI_packet_info mpi_packet_info_ghost;
  

  public:
  };

} // writer

CAVIAR_NAMESPACE_CLOSE

#endif
