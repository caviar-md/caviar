
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

    caviar::Atom_data *atom_data;

    // used in msd calculations
    Domain *domain;

    void dump_energy(int64_t);         // dump energy to file
    void dump_energy(int64_t, double); // dump energy to file

    void dump_temperature(int64_t);         // dump temperature to file
    void dump_temperature(int64_t, double); // dump temperature to file

    void dump_xyz(int64_t);         // dump positions to file in xyz format
    void dump_xyz(int64_t, double); // dump positions to file in xyz format

    void dump_povray(int64_t);         // dump positions to snapshot files in povray format
    void dump_povray(int64_t, double); // dump positions to snapshot files in povray format

    void dump_msd(int64_t);         //
    void dump_msd(int64_t, double); //

    int64_t energy_step, temperature_step, xyz_step, povray_step, msd_step; // number of steps to output data

    int msd_type; // type of atom that should be used in msd calculations.
    int msd_initial_step;

    std::ofstream ofs_energy, ofs_temperature, ofs_xyz, ofs_velocities, ofs_povray; // output files
    std::ofstream ofs_msd;                                                          // mean square distance

    // if true, outputs would be created
    bool output_energy, output_temperature, output_xyz, output_povray, output_msd;

    // dump velocity and acceleration alongside position
    bool output_velocity, output_acceleration;

    std::vector<caviar::Vector<double>> msd_initial_position;

    // records previous wallTime of XYZ dump.
    double wallTimeXyzDump1;

    /////////

  public:
  };

} // writer

CAVIAR_NAMESPACE_CLOSE

#endif
