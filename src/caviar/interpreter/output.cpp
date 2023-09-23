
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

#include "caviar/interpreter/output.h"
#include "caviar/interpreter/communicator.h"
#include "caviar/utility/interpreter_io_headers.h"
#if defined(CAVIAR_WITH_MPI)
#include<mpi.h>
#endif
CAVIAR_NAMESPACE_OPEN
namespace interpreter
{
  Output::Output(CAVIAR *fptr) : Pointers{fptr}
  {
    for (int i = 0; i < 5; ++i)
    {
      output_info[i] = true;
      output_warning[i] = true;
    }
  }

  void Output::info(const char *str, int level, bool endline)
  {
    info(static_cast<std::string>(str), level, endline);
  }
  void Output::info(const std::string &str, int level, bool endline)
  {
    if (output_info[level])
    {
#if defined(CAVIAR_WITH_MPI)
      MPI_Barrier(MPI_COMM_WORLD);
      int me = comm->me;
      if (me == 0)
      {
#endif
        std::cout << "[INF] ";
        std::cout << str;
        if (endline)
          std::cout << std::endl;
        else
          std::cout << std::flush;
#if defined(CAVIAR_WITH_MPI)
      }
#endif
    }
  }

  void Output::info_create(const char *str, int level, bool endline)
  {
    info_create(static_cast<std::string>(str), level, endline);
  }
  void Output::info_create(const std::string &str, int level, bool endline)
  {
    std::string s = "(Create) " + str ;
    info(s, level, endline);
  }

  void Output::info_read(const char *str, int level, bool endline)
  {
    info_read(static_cast<std::string>(str), level, endline);
  }
  void Output::info_read(const std::string &str, int level, bool endline)
  {
    std::string s = "(Call) "+ str + ".read()";
    info(s, level, endline);
  }

  void Output::warning(const char *str, int level, bool endline)
  {
    warning(static_cast<std::string>(str), level, endline);
  }
  void Output::warning(const std::string &str, int level, bool endline)
  {
    if (output_warning[level])
    {
#if defined(CAVIAR_WITH_MPI)
      MPI_Barrier(MPI_COMM_WORLD);
      int me = comm->me;
      if (me == 0)
      {
#endif
        std::cout << "[WRN] ";
        std::cout << str;
        if (endline)
          std::cout << std::endl;
        else
          std::cout << std::flush;
#if defined(CAVIAR_WITH_MPI)
      }
#endif
    }
  }

  void Output::comment(const char *str, bool endline)
  {
    comment(static_cast<std::string>(str), endline);
  }
  void Output::comment(const std::string &str, bool endline)
  {
#if defined(CAVIAR_WITH_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    int me = comm->me;
    if (me == 0)
    {
#endif
      std::cout << str;
      if (endline)
        std::cout << std::endl;
      else
        std::cout << std::flush;
#if defined(CAVIAR_WITH_MPI)
    }
#endif
  }

  bool Output::read(caviar::interpreter::Parser *parser)
  {
    info("output read");
    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      if (string_cmp(t, "info"))
      {
        int level = -1;
        GET_OR_CHOOSE_A_INT(level, "", "")
        bool status = 0;
        GET_A_BOOL(status, "", "")
        if (level > 5 || level < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "'info' level is defined between 0 to 4. To choose all, use 5.");
        if (level == 5)
          for (int i = 0; i < 5; ++i)
            output_info[i] = status;
        else
          output_info[level] = status;
      }
      else if (string_cmp(t, "warning"))
      {
        int level = -1;
        GET_OR_CHOOSE_A_INT(level, "", "")
        bool status = 0;
        GET_A_BOOL(status, "", "")
        if (level > 5 || level < 0)
          error->all(FC_FILE_LINE_FUNC_PARSE, "'warning' level is defined between 0 to 4. To choose all, use 5.");
        if (level == 5)
          for (int i = 0; i < 5; ++i)
            output_warning[i] = status;
        else
          output_warning[level] = status;
      }
      else
        error->all(FC_FILE_LINE_FUNC_PARSE, "Invalid syntax: This output command doesn't exist.");
    }
    return in_file;

    /*
      auto command = parser->get_identifier();
      if (command=="energy") {
        auto steps = parser->get_literal_int();
        energy_step = steps;
        output_energy = true;
      } else if (command=="xyz") {
        auto steps = parser->get_literal_int();
        xyz_step = steps;
        output_xyz = true;
      } else if (command=="povray") {
        auto steps = parser->get_literal_int();
        povray_step = steps;
        output_povray = true;
      } else {
          error->all (FC_FILE_LINE_FUNC_PARSE, "Invalid syntax: This output command doesn't exist.");
      }

      if (! parser->end_of_line ()) error->all (FC_FILE_LINE_FUNC_PARSE, "Invalid syntax"); */
  }

} // interpreter
CAVIAR_NAMESPACE_CLOSE
