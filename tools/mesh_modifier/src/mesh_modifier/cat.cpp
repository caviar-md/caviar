
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

#include "mesh_modifier.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <iterator>
#include <iomanip>

namespace mesh_modifier
{

  void Mesh_modifier::cat_by_string(const std::string &filename)
  {
    std::ifstream ifs;
    ifs.open(filename.c_str());
    while (!ifs.eof())
    {
      // char c;
      std::string c;
      ifs >> c;
      std::cout << c << " ";
      if (c == "-1")
        std::cout << "\n";
    }
    ifs.close();
  }

  void Mesh_modifier::cat_by_line(const std::string &filename)
  {
    std::ifstream ifs;
    ifs.open(filename.c_str());
    std::string line;
    while (std::getline(ifs, line))
    {
      std::cout << line << std::endl;
    }
    ifs.close();
  }

  void Mesh_modifier::cat_by_line_into_token(const std::string &filename)
  {
    std::ifstream ifs;
    ifs.open(filename.c_str());
    std::string line;
    while (std::getline(ifs, line))
    {
      // std::cout << line << std::endl;
      std::vector<std::string> tokens;
      std::istringstream iss(line);
      std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(tokens));
      for (unsigned int i = 0; i < tokens.size(); ++i)
        std::cout << tokens[i] << " ";
      std::cout << "\n";
    }
    ifs.close();
  }

  void Mesh_modifier::cat(const std::string &filename, bool unsupported)
  {
    std::ifstream ifs;
    ifs.open(filename.c_str());
    bool ignore_section = false;
    bool in_section = false;
    int udn_code = 0;
    while (!ifs.eof())
    {
      if (!in_section)
      {
        std::string c;
        ifs >> c;
        if (c == "#")
          ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if (c == "-1")
        {
          std::cout << std::right << std::setw(6) << "-1"
                    << "\n";
          in_section = true;
          ifs >> udn_code;
          std::cout << std::right << std::setw(6) << udn_code << "\n";
        }
      }
      else
      {
        if (udn_code == 2411)
          cat_udn_2411(ifs);
        else if (udn_code == 2412)
          cat_udn_2412(ifs);
        else if (udn_code == 2467)
          cat_udn_2467(ifs);
        else
        {

          if (unsupported)
            cat_udn_unsupported(ifs, udn_code);
          else
            cat_udn_ignore(ifs, udn_code);
        }
        in_section = false;
      }
    }
    ifs.close();
  }

  void Mesh_modifier::cat_udn_ignore(std::ifstream &ifs, int udn_code)
  {
    std::cout << "Warning: Unsupported udn " << udn_code
              << ". Ignoring the section.\n";
    while (true)
    {
      std::string c;
      ifs >> c;
      if (c == "-1")
        return;
    }
  }

  void Mesh_modifier::cat_udn_unsupported(std::ifstream &ifs, int udn_code)
  {
    std::cout << "Warning: Unsupported udn " << udn_code << "\n";
    // std::cout << std::right << std::setw(6) << udn_code << "\n";
    while (true)
    {
      std::string line;
      std::getline(ifs, line);
      std::vector<std::string> tokens;
      std::istringstream iss(line);
      std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(tokens));
      if (tokens.size() > 0)
      {
        if (tokens[0] == "-1")
        {
          std::cout << std::right << std::setw(6) << "-1"
                    << "\n";
          return;
        }
        std::cout << line << "\n";
      }
    }
  }

  void Mesh_modifier::cat_udn_2411(std::ifstream &ifs)
  {

    while (true)
    {
      int tmp;
      ifs >> tmp;
      if (tmp == -1)
      {
        std::cout << std::right << std::setw(6) << "-1"
                  << "\n";
        return;
      }

      int dummy[2];
      ifs >> dummy[0] >> dummy[1] >> dummy[2];
      //      std::cout << tmp << " " <<  dummy[0] << " " << dummy[1] << " " << dummy[2] << "\n";
      std::cout << std::right << std::setw(10) << tmp
                << std::setw(10) << dummy[0]
                << std::setw(10) << dummy[1]
                << std::setw(10) << dummy[2] << "\n";

      double pos[2];
      ifs >> pos[0] >> pos[1] >> pos[2];
      std::cout << std::scientific << std::setprecision(16) << std::uppercase
                << std::setw(25) << pos[0]
                << std::setw(25) << pos[1]
                << std::setw(25) << pos[2] << "\n";
      std::cout << std::fixed;
    }
  }

  void Mesh_modifier::cat_udn_2412(std::ifstream &ifs)
  {

    while (true)
    {
      int tmp;
      ifs >> tmp;
      if (tmp == -1)
      {
        std::cout << std::right << std::setw(6) << "-1"
                  << "\n";
        return;
      }

      int dummy[5];
      ifs >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> dummy[4];

      // std::cout << tmp << " " <<  dummy[0] << " " << dummy[1] << " " << dummy[2]
      //           << " " <<  dummy[3] << " " << dummy[4] << "\n";

      std::cout << std::right << std::setw(10) << tmp
                << std::setw(10) << dummy[0]
                << std::setw(10) << dummy[1]
                << std::setw(10) << dummy[2]
                << std::setw(10) << dummy[3]
                << std::setw(10) << dummy[4] << "\n";

      int FE_Id = dummy[0];
      bool beam_type = (FE_Id == 11 || FE_Id == 21 || FE_Id == 22 || FE_Id == 23 || FE_Id == 24);
      unsigned num_of_elements = dummy[4];

      if (beam_type)
      { // beam elements

        int dummy_[2];
        ifs >> dummy_[0] >> dummy_[1] >> dummy_[2];
        // std::cout <<  dummy_[0] << " " << dummy_[1] << " " << dummy_[2]  << "\n";
        std::cout << std::setw(10) << dummy_[0]
                  << std::setw(10) << dummy_[1]
                  << std::setw(10) << dummy_[2] << "\n";

        int field[num_of_elements];
        for (unsigned int i = 0; i < num_of_elements; ++i)
          ifs >> field[i];
        for (unsigned int i = 0; i < num_of_elements; ++i)
        {
          // std::cout << field[i] << " ";
          std::cout << std::setw(10) << field[i];
        }
        std::cout << "\n";
      }
      else
      {

        int field[num_of_elements];
        for (unsigned int i = 0; i < num_of_elements; ++i)
          ifs >> field[i];
        for (unsigned int i = 0; i < num_of_elements; ++i)
        {
          // std::cout << field[i] << " ";
          std::cout << std::setw(10) << field[i];
        }
        std::cout << "\n";
      }
    }
  }

  void Mesh_modifier::cat_udn_2467(std::ifstream &ifs)
  {

    while (true)
    {

      int tmp;
      ifs >> tmp;
      if (tmp == -1)
        return;

      int dummy[7];
      ifs >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> dummy[4] >> dummy[5] >> dummy[6];

      //      std::cout << tmp << " " <<  dummy[0] << " " << dummy[1] << " "
      //                << dummy[2]   << " " << dummy[3] << " " << dummy[4] << " "
      //                << dummy[5]   << " " << dummy[6] << "\n";

      std::cout << std::right << std::setw(10) << tmp
                << std::setw(10) << dummy[0]
                << std::setw(10) << dummy[1]
                << std::setw(10) << dummy[2]
                << std::setw(10) << dummy[3]
                << std::setw(10) << dummy[4]
                << std::setw(10) << dummy[5]
                << std::setw(10) << dummy[6] << "\n";

      int group_name;
      ifs >> group_name;
      std::cout << group_name << "\n";

      bool print_eol = false;
      for (int i = 0; i < dummy[6]; ++i)
      {
        int dummy[4];
        ifs >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3];

        //        std::cout << dummy[0] << " " << dummy[1] << " " << dummy[2] << " "
        //                  << dummy[3] ;

        std::cout << std::right
                  << std::setw(10) << dummy[0]
                  << std::setw(10) << dummy[1]
                  << std::setw(10) << dummy[2]
                  << std::setw(10) << dummy[3];

        if (print_eol)
        {
          std::cout << "\n";
          print_eol = false;
        }
        else
        {
          //          std::cout << " ";
          print_eol = true;
        }
      }
      if (print_eol)
        std::cout << "\n";
    }
  }
}
