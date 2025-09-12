
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

#ifndef CAVIAR_INTERPRETER_ERROR_H
#define CAVIAR_INTERPRETER_ERROR_H

#include "caviar/utility/pointers.h"

#define FC_FILE_LINE_FUNC_LINE_COL __FILE__, __LINE__, __func__, line, col
#define FC_FILE_LINE_FUNC_PARSE __FILE__, __LINE__, __func__, parser->line, parser->col
#define FC_FILE_LINE_FUNC __FILE__, __LINE__, __func__

CAVIAR_NAMESPACE_OPEN
namespace interpreter
{

  /**
   * This class does all the error handling outputs.
   * It has different types of error calls.
   *
   */
  class Error : public Pointers
  {
  public:
    /**
     * Constructor.
     */
    Error(CAVIAR *);

    void all(const std::string &);

    void all(const char *, int, const char *, const std::string &, unsigned int, const char *);
    void one(const char *, int, const char *, const std::string &, unsigned int, const char *);

    void all(const char *, int, const char *, const std::string &, unsigned int, const std::string &);
    void one(const char *, int, const char *, const std::string &, unsigned int, const std::string &);

    void all(const char *, int, const char *, const std::string &);
    void one(const char *, int, const char *, const std::string &);
  };
} // interpreter
CAVIAR_NAMESPACE_CLOSE

#endif
