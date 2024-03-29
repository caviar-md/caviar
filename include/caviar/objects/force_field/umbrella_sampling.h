
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_UMBRELLA_SAMPLING_H
#define CAVIAR_OBJECTS_FORCEFIELD_UMBRELLA_SAMPLING_H

#include "caviar/objects/force_field.h"

CAVIAR_NAMESPACE_OPEN

namespace force_field {

/**
 * This class does a spring force-field on the molecular bonds
 *  
 */
class Umbrella_sampling : public Force_field {
public:
  Umbrella_sampling (class CAVIAR *);
  ~Umbrella_sampling () {};

  //double energy();

  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:
  void init_production();

  void production_function();

  void finish_production();

  void init_metadata();

  void finish_metadata();

  /**
   * If true, the metadata is written in the file
  */
  bool metadata_mode = false;

  /**
   * If true, the data is written in the file
  */
  bool production_mode = false;

  /**
   * The elastic coefficient of the force
  */
  double elastic_coef = 0;

  /**
   * The atom_id in which the force is acted on
  */
  int atom_id = -1;

  /**
   * The anchor point position, i.e., the constraint position
  */
  caviar::Vector<double> position {0,0,0};

  /**
   * The spring 
  */
  caviar::Vector<double> dr;

  /**
   * Count number of files which is produced
  */
  int file_counter = 0;

  /**
   * Count the number of steps locally in order to use in skipping frames
  */
  int step_counter = 0;

  /**
   * Count the number of steps locally in order to use in skipping frames
  */
  char reaction_coordinate = 'x';

  /**
   * Every number of 'step' the data is written in the file
  */
  int step = 1;

  /**
   * set temperature to be outputted in the file
  */
  double temperature = 1;

  /**
   * set temperature 
  */
  bool temperature_is_set = false;

  /**
   * Make an absolute path for the metadata file
  */
  void fix_file_prefix();

  double position_min = 0;
  double position_mean = 0;
  int position_mean_count = 0;
  double position_max = 0;

  std::vector<double> session_data;

  std::ofstream ofs_data;
  std::ofstream ofs_metadata;
  std::ofstream ofs_stat;

  std::string file_name_data = "o_umbr_";
  std::string file_name_metadata = "o_umbr_metadata";
  std::string file_name_stat = "o_umbr_stat";
  std::string file_prefix = "";

};

} //force_field

CAVIAR_NAMESPACE_CLOSE

#endif
