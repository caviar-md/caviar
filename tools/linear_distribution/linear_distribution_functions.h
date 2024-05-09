
//========================================================================
//
// Copyright (C) 2024 by Morad Biagooi and Ehsan Nedaaee Oskoee.
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

#ifndef LINEAR_DISTRIBUTOIN_FUNCTION
#define LINEAR_DISTRIBUTOIN_FUNCTION

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <sys/stat.h> // For mkdir

// Function to check if a directory exists
inline bool directoryExists(const std::string &dirName)
{
    struct stat info;
    if (stat(dirName.c_str(), &info) != 0)
        return false;
    else if (info.st_mode & S_IFDIR) // S_ISDIR() doesn't exist on my windows
        return true;
    else
        return false;
}

// Function to generate a unique output folder name
inline std::string generateOutputFolderName(const std::string &xyz_filename)
{
    std::string folderName = "out_" + xyz_filename;
    int counter = 1;
    while (directoryExists(folderName))
    {
        std::stringstream ss;
        ss << "out_" << xyz_filename << "_" << counter++;
        folderName = ss.str();
    }
    return folderName;
}

struct Lim
{
    bool enable;
    double min, max;
};

struct Lim3D
{
    Lim x,y,z;
};

struct Vec3D
{
    double x, y, z;

    int Normalize()
    {
        double norm = std::sqrt(x * x + y * y + z * z);
        if (norm == 0)
        {
            std::cout << "Error: Cannot normalize a vector of length 0" << std::endl;
            return 1;
        }
        x = x / norm;
        y = y / norm;
        z = z / norm;
        return 0;
    }
};

struct Atom
{
    Vec3D pos;
    int type;
};

struct Snapshot
{
    std::vector<Atom> atoms;
};

std::vector<Snapshot> readXYZFile(const std::string &filename);

std::vector<double> calculateAverageDistribution(const std::vector<Snapshot> &snapshots,
                                                 int atomType,
                                                 Vec3D start,
                                                 Vec3D direction,
                                                 int numBins, double totalLength);

std::vector<double> calculateLinearDistribution(const Snapshot &snapshot,
                                                int atomType,
                                                Vec3D start,
                                                Vec3D direction,
                                                int numBins, double totalLength);

#endif