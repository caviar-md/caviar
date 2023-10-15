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
//
// some part of this file (file_exists_ functions) are developed by editing the best answer to the following questions
// https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exists-using-standard-c-c11-14-17-c
//

#include <sys/stat.h>
// #include <unistd.h>
#include <string>
// #include <fstream>

// inline bool file_exists_0(const std::string &name)
// {
//     std::ifstream f(name.c_str());
//     return f.good();
// }

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

// inline bool file_exists_2(const std::string &name)
// {
//     return (access(name.c_str(), F_OK) != -1);
// }

// inline bool file_exists_3(const std::string &name)
// {
//     struct stat buffer;
//     return (stat(name.c_str(), &buffer) == 0);
// }

inline std::string directory_of_file(const std::string &file_path)
{
    // Find the last occurrence of '/' or '\' to determine the directory
    size_t found = file_path.find_last_of("/\\");

    if (found != std::string::npos)
    {
        std::string directory = file_path.substr(0, found);
        return directory;
    }
    else
    {
        return "";
    }
}

inline std::string join_path(const std::string &directory, const std::string &fileName)
{
    if (directory == "")
        return fileName;

    if (!directory.empty() && directory.back() != '/' && directory.back() != '\\')
    {
        // Add a path separator if it's missing
        return directory + '/' + fileName;
    }
    else
    {
        return directory + fileName;
    }
}

#include <unistd.h>
#include <limits.h>

inline std::string get_current_directory()
{
   char cwd[PATH_MAX];
   std::string st = "";
   if (getcwd(cwd, sizeof(cwd)) != NULL) 
   {
       //printf("Current working dir: %s\n", cwd);
       st = cwd;
   } 
   else 
   {
       //perror("getcwd() error");
       //return 1;
   }
   return st;
}

inline bool is_absolute_path(const std::string &path)
{
    return !path.empty() && path[0] == '/';
}