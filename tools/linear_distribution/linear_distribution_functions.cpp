
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

#include "linear_distribution_functions.h"

extern Lim3D lim3d;

// Function to read XYZ file
std::vector<Snapshot> readXYZFile(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<Snapshot> snapshots;
    Snapshot snapshot;
    std::string line;
    int numAtoms;
    bool first_line = true;
    int atomCount = 0;
    while (std::getline(file, line))
    {
        // std::cout << line << std::endl;

        std::istringstream iss(line);

        if (first_line)
        {
            atomCount = 0;
            snapshot.atoms.clear();

            first_line = false;
            if (iss >> numAtoms)
            {
                continue;
            }
            else
            {
                break;
            }
        }

        if (line.find("Atom") != std::string::npos)
        {
            continue; // Ignore lines containing ignoreText
        }

        Atom atom;
        if (!(iss >> atom.type >> atom.pos.x >> atom.pos.y >> atom.pos.z))
        {
            std::cerr << "Error: Invalid format in XYZ file." << std::endl;
            exit(EXIT_FAILURE);
        }
        atomCount++;
        // std::cout << "numAtoms: " << numAtoms << " atomCount: " << atomCount << std::endl;
        snapshot.atoms.push_back(atom);

        if (atomCount == numAtoms)
        {
            snapshots.push_back(snapshot);
            first_line = true;
        }
    }

    file.close();
    return snapshots;
}

std::vector<double> calculateLinearDistribution(const Snapshot &snapshot,
                                                int atomType,
                                                Vec3D start,
                                                Vec3D direction,
                                                int numBins, double totalLength)
{
    double binLength = totalLength / numBins;

    std::vector<double> distribution(numBins, 0.0);
    int count = 0;

    for (const auto &atom : snapshot.atoms)
    {
        if (atom.type != atomType)
            continue;

        if (lim3d.x.enable)
        {
            if (atom.pos.x < lim3d.x.min || atom.pos.x > lim3d.x.max)
                continue;
        }
        if (lim3d.y.enable)
        {
            if (atom.pos.y < lim3d.y.min || atom.pos.y > lim3d.y.max)
                continue;
        }
        if (lim3d.z.enable)
        {
            if (atom.pos.z < lim3d.z.min || atom.pos.z > lim3d.z.max)
                continue;
        }

        double dx = atom.pos.x - start.x;
        double dy = atom.pos.y - start.y;
        double dz = atom.pos.z - start.z;
        double projection = dx * direction.x + dy * direction.y + dz * direction.z;
        if (projection >= 0 && projection <= totalLength)
        {
            int bin = static_cast<int>(projection / binLength);
            distribution[bin]++;
            count++;
        }
    }

    if (count > 0)
    {
        // Normalize distribution
        for (double &value : distribution)
        {
            value /= count;
        }
    }

    return distribution;
}

std::vector<double> calculateAverageDistribution(const std::vector<Snapshot> &snapshots,
                                                 int atomType,
                                                 Vec3D start,
                                                 Vec3D direction,
                                                 int numBins, double totalLength)
{

    double binLength = totalLength / numBins;

    std::vector<double> averageDistribution(numBins, 0.0);

    for (const auto &snapshot : snapshots)
    {
        std::vector<double> distribution = calculateLinearDistribution(snapshot, atomType, start, direction, numBins, totalLength);
        for (int i = 0; i < numBins; ++i)
        {
            averageDistribution[i] += distribution[i];
        }
    }

    // Normalize average distribution
    int numSnapshots = snapshots.size();
    for (double &value : averageDistribution)
    {
        value /= numSnapshots;
    }

    return averageDistribution;
}
