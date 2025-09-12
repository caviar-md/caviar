
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
#include "caviar/CAVIAR.hpp"
#include "caviar/interpreter/all.hpp"
#include "caviar/interpreter/error.hpp"
#if defined(CAVIAR_WITH_MPI)
#include <mpi.h>
#endif
namespace caviar
{
  namespace interpreter
  {
    Error::Error(CAVIAR *fptr) : caviar_{fptr},
                                   comm{fptr->comm},
                                   error{fptr->error},
                                   output{fptr->output},
                                   input{fptr->input},
                                   object_handler{fptr->object_handler},
                                   object_container{fptr->object_container},
                                   object_creator{fptr->object_creator},
                                   log{fptr->log},
                                   in{fptr->in},
                                   out{fptr->out},
                                   err{fptr->err},
                                   log_flag{fptr->log_flag},
                                   out_flag{fptr->out_flag},
                                   err_flag{fptr->err_flag} {}
Error::~Error(){}
    // All procs must call this else there would be a deadlock
    void Error::all(const std::string &str)
    {
      err << "[ERR] " << str << std::endl;
      exit(1);
    }
    void Error::verify_settings() {}
    void Error::all(const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const char *str)
    {
      int me = 0;
#ifdef CAVIAR_WITH_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &me);
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      if (me == 0)
      {
        if (err_flag)
        {
          err << "[ERR] " << str << std::endl;
          err << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
          err << parsing_line << std::endl;
          for (unsigned int i = 0; i < col; ++i)
            err << ' ';
          err << '^' << std::endl;
        }
        if (log_flag)
        {
          log << "[ERR] " << str << std::endl;
          log << " '" << file << ':' << line << " in '" << func << "'." << std::endl;
          log << parsing_line << std::endl;
          for (unsigned int i = 0; i < col; ++i)
            log << ' ';
          log << '^' << std::endl;
        }
        if (log_flag)
          log.close();
      }
#ifdef CAVIAR_WITH_MPI
      MPI_Finalize();
#endif
      exit(1);
    }

    // One proc calling this will abort all

    void Error::one(const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const char *str)
    {
#ifdef CAVIAR_WITH_MPI
      int me;
      MPI_Comm_rank(MPI_COMM_WORLD, &me);
      if (err_flag)
      {
        err << "[ERR] " << str << std::endl;
        err << " MPI rank " << me << std::endl;
        err << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;

        err << parsing_line << std::endl;
        for (unsigned i = 0; i < col; ++i)
          err << ' ';
        err << '^' << std::endl;
      }
      if (log_flag)
      {
        log << "[ERR] " << str << std::endl;
        err << " MPI rank " << me << std::endl;
        log << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
        log << parsing_line << std::endl;
        for (unsigned i = 0; i < col; ++i)
          log << ' ';
        log << '^' << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
#else
      all(file, line, func, parsing_line, col, str);
#endif
    }

    void Error::all(const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const std::string &str)
    {
      int me = 0;
#ifdef CAVIAR_WITH_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &me);
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      //  /*
      if (me == 0)
      {
        if (err_flag)
        {
          err << "[ERR] " << str << std::endl;
          err << " MPI rank " << me << std::endl;
          err << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
          err << parsing_line << std::endl;
          for (unsigned int i = 0; i < col; ++i)
            err << ' ';
          err << '^' << std::endl;
        }
        if (log_flag)
        {
          log << "[ERR] " << str << std::endl;
          err << " MPI rank " << me << std::endl;
          log << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
          log << parsing_line << std::endl;
          for (unsigned int i = 0; i < col; ++i)
            log << ' ';
          log << '^' << std::endl;
        }
        if (log_flag)
          log.close();
      }
//  */
#ifdef CAVIAR_WITH_MPI
      MPI_Finalize();
#endif
      exit(1);
    }

    void Error::one(const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const std::string &str)
    {
#ifdef CAVIAR_WITH_MPI
      int me;
      MPI_Comm_rank(MPI_COMM_WORLD, &me);
      if (err_flag)
      {
        err << "[ERR] " << str << std::endl;
        err << " MPI rank " << me << std::endl;
        err << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
        err << parsing_line << std::endl;
        for (unsigned i = 0; i < col; ++i)
          err << ' ';
        err << '^' << std::endl;
      }
      if (log_flag)
      {
        log << "[ERR] " << str << std::endl;
        err << " MPI rank " << me << std::endl;
        log << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
        log << parsing_line << std::endl;
        for (unsigned i = 0; i < col; ++i)
          log << ' ';
        log << '^' << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
#else
      all(file, line, func, parsing_line, col, str);
#endif
    }

    void Error::all(const char *file, int line, const char *func, const std::string &str)
    {
      int me = 0;
#ifdef CAVIAR_WITH_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &me);
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      if (me == 0)
      {
        if (err_flag)
        {
          err << "[ERR] " << str << std::endl;
          err << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
        }
        if (log_flag)
        {
          log << "[ERR] " << str << std::endl;
          log << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
        }
        if (log_flag)
          log.close();
      }

#ifdef CAVIAR_WITH_MPI
      MPI_Finalize();
#endif
      exit(1);
    }

    void Error::one(const char *file, int line, const char *func, const std::string &str)
    {
#ifdef CAVIAR_WITH_MPI
      int me;
      MPI_Comm_rank(MPI_COMM_WORLD, &me);
      if (err_flag)
      {
        err << "[ERR] " << str << std::endl;
        err << " MPI rank " << me << std::endl;
        err << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
      }
      if (log_flag)
      {
        log << "[ERR] " << str << std::endl;
        log << " MPI rank " << me << std::endl;
        log << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
#else
      all(file, line, func, str);
#endif
    }
  } // interpreter
}
