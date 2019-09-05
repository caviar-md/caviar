
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

#include "caviar/interpreter/error.h"

namespace caviar {
namespace interpreter {
Error::Error (CAVIAR *fptr) : Pointers{fptr} {}


// All procs must call this else there would be a deadlock
void Error::all (const std::string &str) {
  err << "  - ERROR occured!"  <<std::endl;
  err << "  - ERROR massage: " << str << std::endl;
  exit(1);
}

void Error::all (const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const char *str) {
  int me = 0;
#ifdef CAVIAR_WITH_MPI
  MPI_Comm_rank (mpi_comm, &me);
  MPI_Barrier (mpi_comm);
#endif
  std::cout << "  - ERROR occured!"  <<std::endl;
  if (me == 0) {
    if (err_flag) {
      err << "  - ERROR massage: " << str << std::endl;
      err << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;
      err << parsing_line << std::endl;
      for (unsigned int i=0; i<col; ++i) err << ' ';
      err << '^' << std::endl;
    }
    if (log_flag) {
      log << "  - ERROR massage: " << str << std::endl;
      log << "  - ERROR call: at " << file << ':' << line << " in '" << func << "'." << std::endl;
      log << parsing_line << std::endl;
      for (unsigned int i=0; i<col; ++i) log << ' ';
      log << '^' << std::endl;
    }
    if (log_flag) log.close(); 
  }
#ifdef CAVIAR_WITH_MPI
  MPI_Finalize ();
#endif
  exit (1);
}

// One proc calling this will abort all

void Error::one (const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const char *str) {
#ifdef CAVIAR_WITH_MPI
  int me;
  MPI_Comm_rank (mpi_comm, &me);
  std::cout << "  - ERROR occured!"  <<std::endl;
  if (err_flag) {
    err << "  - ERROR on proc " << me << std::endl;
    err << "  - ERROR massage: " << str << std::endl;
    err << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;
    err << parsing_line << std::endl;
    for (unsigned i=0; i<col; ++i) err << ' ';
    err << '^' << std::endl;
  }
  if (log_flag) {
    log << "  - ERROR on proc " << me << std::endl;
    log << "  - ERROR massage: " << str << std::endl;
    log << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;
    log << parsing_line << std::endl;
    for (unsigned i=0; i<col; ++i) log << ' ';
    log << '^' << std::endl;
  }
  MPI_Abort (mpi_comm, 1);
#else
  all (file, line, func, parsing_line, col, str);
#endif
}

void Error::all (const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const std::string &str) {
  int me = 0;
#ifdef CAVIAR_WITH_MPI
  MPI_Comm_rank (mpi_comm, &me);
  MPI_Barrier (mpi_comm);
#endif
  std::cout << "ERROR occured!"  <<std::endl;
//  /*
  if (me == 0) {
    if (err_flag) {
      err << "  - ERROR massage: " << str << std::endl;
      err << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;
      err << parsing_line << std::endl;
      for (unsigned int i=0; i<col; ++i) err << ' ';
      err << '^' << std::endl;
    }
    if (log_flag) {
      log << "  - ERROR massage: " << str << std::endl;
      log << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;
      log << parsing_line << std::endl;
      for (unsigned int i=0; i<col; ++i) log << ' ';
      log << '^' << std::endl;
    }
    if (log_flag) log.close(); 
  }
//  */
#ifdef CAVIAR_WITH_MPI
  MPI_Finalize ();
#endif
  exit (1);
}

void Error::one (const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const std::string &str) {
#ifdef CAVIAR_WITH_MPI
  int me;
  MPI_Comm_rank (mpi_comm, &me);
  std::cout << "  - ERROR occured!"  <<std::endl;
  if (err_flag) {
    err << "  - ERROR on proc " << me << std::endl;
    err << "  - ERROR massage: " << str << std::endl;
    err << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;
    err << parsing_line << std::endl;
    for (unsigned i=0; i<col; ++i) err << ' ';
    err << '^' << std::endl;
  }
  if (log_flag) {
    log << "  - ERROR on proc " << me << std::endl;
    log << "  - ERROR massage: " << str << std::endl;
    log << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;
    log << parsing_line << std::endl;
    for (unsigned i=0; i<col; ++i) log << ' ';
    log << '^' << std::endl;
  }
  MPI_Abort (mpi_comm, 1);
#else
  all (file, line, func, parsing_line, col, str);
#endif
}

void Error::all (const char *file, int line, const char *func, const std::string &str) {
  int me = 0;
#ifdef CAVIAR_WITH_MPI
  MPI_Comm_rank (mpi_comm, &me);
  MPI_Barrier (mpi_comm);
#endif
  std::cout << "  - ERROR occured!"  <<std::endl;

  if (me == 0) {
    if (err_flag) {
      err << "  - ERROR massage: " << str << std::endl;
      err << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;
    }
    if (log_flag) {
      log << "  - ERROR massage: " << str << std::endl;
      log << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;
    }
    if (log_flag) log.close(); 
  }

#ifdef CAVIAR_WITH_MPI
  MPI_Finalize ();
#endif
  exit (1);
}

void Error::one (const char *file, int line, const char *func, const std::string &str) {
#ifdef CAVIAR_WITH_MPI
  int me;
  MPI_Comm_rank (mpi_comm, &me);
  std::cout << "  - ERROR occured!"  <<std::endl;
  if (err_flag) {
    err << "  - ERROR on proc " << me << std::endl;
    err << "  - ERROR massage: " << str << std::endl;
    err << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;

  }
  if (log_flag) {
    log << "  - ERROR on proc " << me << std::endl;
    log << "  - ERROR massage: " << str << std::endl;
    log << "  - ERROR call: at '" << file << ':' << line << "' in '" << func << "'." << std::endl;

  }
  MPI_Abort (mpi_comm, 1);
#else
  all (file, line, func, str);
#endif
}
} //interpreter
} // namespace caviar

