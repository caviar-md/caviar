
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

#include "caviar/objects/postprocess/potential_sampler.h"
#include "caviar/utility/interpreter_io_headers.h"
#include "caviar/objects/atom_data.h"
#include "caviar/objects/force_field.h"
#include "caviar/objects/md_simulator.h"
#include "caviar/objects/unique/grid_1d.h"

CAVIAR_NAMESPACE_OPEN

namespace postprocess
{

  Potential_sampler::Potential_sampler(CAVIAR *fptr) : Postprocess{fptr}
  {
    FC_OBJECT_INITIALIZE_INFO
    step_start = 0;
    step_end = -1;
    output_file_name = "o_potential_sampler";
    read_velocity = false;
    atom_data = nullptr;
    md_simulator = nullptr;
    grid_x = nullptr;
    grid_y = nullptr;
    grid_z = nullptr;
  }

  Potential_sampler::~Potential_sampler() {}

  bool Potential_sampler::read(caviar::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO

    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR

      if (t == "run")
      {
        run();
        break;
      }
      else if (string_cmp(t, "input_xyz_file_name"))
      {
        const auto token = parser->get_val_token();
        const auto file_name = token.string_value;
        input_xyz_file_name = file_name;
      }
      else if (string_cmp(t, "output_file_name"))
      {
        const auto token = parser->get_val_token();
        const auto file_name = token.string_value;
        output_file_name = file_name;
      }
      else if (string_cmp(t, "add_force_field") || string_cmp(t, "force_field"))
      {
        FIND_OBJECT_BY_NAME(force_field, it)
        force_field.push_back(object_container->force_field[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_x"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_x = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_y"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_y = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_z"))
      {
        FIND_OBJECT_BY_NAME(unique, it)
        grid_z = static_cast<unique::Grid_1D *>(object_container->unique[it->second.index]);
      }
      else if (string_cmp(t, "set_md_simulator") || string_cmp(t, "md_simulator"))
      {
        FIND_OBJECT_BY_NAME(md_simulator, it)
        md_simulator = object_container->md_simulator[it->second.index];
      }
      else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else
        ASSIGN_INT(step_start, "GRID_1D Read: ", "")
      else ASSIGN_INT(step_end, "GRID_1D Read: ", "") else ASSIGN_INT(step_increment, "GRID_1D Read: ", "") else error->all(FC_FILE_LINE_FUNC_PARSE, "Random_1D Read: Unknown variable or command ");
    }
    return in_file;
    ;
  }

  void Potential_sampler::run()
  {
    std::cout << "[INF] Potential_sampler::run () " << std::endl;

    FC_NULLPTR_CHECK(atom_data)
    FC_NULLPTR_CHECK(md_simulator)

    set_positions_vectors();

    std::cout << "r1 : input_xyz_file_name : " << input_xyz_file_name << std::endl;
    atom_data->initialize_reading_xyz_frames(input_xyz_file_name);

    step_current = -1;

    ofs_out.open(output_file_name.c_str());

    while (true)
    {
      step_current++;

      if (step_end != -1 && step_current > step_end)
        break;

      bool set_step = (step_current >= step_start);

      int read_result = atom_data->read_next_xyz_frame(set_step, read_velocity);

      if (read_result == -1)
        break;

      if (set_step)
      {
        std::cout << "[INF] Potential_sampler::run: processing step " << step_current << std::endl;

        md_simulator->step(step_current); // check the position of water after applying bond constraints

        sample_potential(); // done
      }
    }

    atom_data->finalize_reading_xyz_frames();

    ofs_out.close();
  }

  int Potential_sampler::read_next_frame(bool, bool)
  {
    // atom_data->read_next_xyz_frame(set_frame, read_velocity);
    return 0;
  }

  void Potential_sampler::sample_potential()
  {

    // int no_forces = force_field.size();

    int i = -1;
    for (auto &pos : sampling_position)
    {
      i++;

      ofs_out << step_current << " " << i << " ";

      ofs_out << sampling_position_index[i].x << " "
              << sampling_position_index[i].y << " "
              << sampling_position_index[i].z << " ";

      ofs_out << pos.x << " "
              << pos.y << " "
              << pos.z << " ";

      double sum_pots = 0.0;

      for (auto &f : force_field)
      {
        double pot = f->potential(pos);

        sum_pots += pot;

        ofs_out << pot << " ";
      }

      ofs_out << sum_pots << "\n";
    }
  }

  void Potential_sampler::set_positions_vectors()
  {

    FC_NULLPTR_CHECK(grid_x)
    FC_NULLPTR_CHECK(grid_y)
    FC_NULLPTR_CHECK(grid_z)

    sampling_position.clear();
    sampling_position_index.clear();

    int num_of_positions = grid_x->no_points() * grid_y->no_points() * grid_z->no_points();

    sampling_position.reserve(num_of_positions);
    sampling_position_index.reserve(num_of_positions);

    grid_x->reset();
    grid_y->reset();
    grid_z->reset();

    for (int i = 0; i < (int)grid_x->no_points(); ++i)
    {
      double x = grid_x->give_point(i);
      for (int j = 0; j < (int)grid_y->no_points(); ++j)
      {
        double y = grid_y->give_point(j);
        for (int k = 0; k < (int)grid_z->no_points(); ++k)
        {
          double z = grid_z->give_point(k);
          sampling_position.push_back(caviar::Vector<double>{x, y, z});
          sampling_position_index.push_back(caviar::Vector<int>{i, j, k});
        }
      }
    }
  }

  /* // -------------------------- grid_1d functions

    output->info("Grid_1d Generate: ");
    if (generated == true)
      error->all("Grid_1D: Generate: cannot be generated twice. ");
    generated = true;
    if (segment>0 && increment>0)
      error->all("Grid_1D: Generate: Assigning both segment and increment is not possible. ");
    if (segment<0 && increment<0)
      error->all("Grid_1D: Generate: Assign one of segment or increment. ");
    if (min > max)
      error->all("Grid_1D: Generate: min has to be smaller than max. ");

    if (segment<0) { // by increment
      by_increment = true;
      segment = int ((max-min)/increment);
    }

    if (increment<0) { // by segment
      by_segment = true;
      increment = (max - min)/double(segment);
    }

  unsigned int Grid_1D::no_points () {
    if (by_segment) {
      if (segment == 0) // segment means that there has to be at least two points
        return 0;
      return segment + 1;

    }  else // by increment. It should have at least one point
      return segment + 1;
  }

  double Grid_1D::give_point () {
    double val = min + no_given_points * increment;
    ++no_given_points;
    if (by_segment) {
      if (no_given_points > segment) return max;
      else return val;
    } else {
      return val;
    }
  }

  double Grid_1D::give_point (int i) {
    double val = min + i * increment;
    if (by_segment) {
      //if (i == segment) return max; // XXX
      //else return val;
      return val;
    } else {
      return val;
    }
  }
  */
} // postprocess

CAVIAR_NAMESPACE_CLOSE
