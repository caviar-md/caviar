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

#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <sstream>

inline bool file_exists_0(const std::string &name)
{
    std::ifstream f(name.c_str());
    return f.good();
}

inline bool file_exists_1(const std::string &name)
{
    if (FILE *file = fopen(name.c_str(), "r"))
    {
        fclose(file);
        return true;
    }
    else
    {
        return false;
    }
}

inline bool file_exists_2(const std::string &name)
{
    return (access(name.c_str(), F_OK) != -1);
}

inline bool file_exists_3(const std::string &name) // best performance
{
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}
