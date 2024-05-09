
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <cmath>

#include "linear_distribution_functions.h"

Lim3D lim3d;

int main(int argc, char *argv[])
{

  // =========================================================
  // Default Input Parameters
  // =========================================================

  std::string xyz_filename = "o_xyz.xyz"; // Path to the XYZ file containing atomic coordinates
  Vec3D start{0.0, 0.0, 0.0};             // Starting position for linear distribution calculation (x, y, z)
  Vec3D direction{1.0, 0.0, 0.0};         // Direction of linear distribution calculation (x, y, z)
  double totalLength = 10.0;              // Total length of the distribution
  int numBins = 50;                       // Number of bins in the length
  int atomTypeStart = 0;                  // Start number of atom type range to be averaged on
  int atomTypeEnd = 4;                    // End number of atom type range to be averaged on

  lim3d.x.enable = false; // Limit the atom positions in the interval in the x-direction.
  lim3d.x.min = 0;
  lim3d.x.max = 0;

  lim3d.y.enable = false; // Limit the atom positions in the interval in the y-direction.
  lim3d.y.min = 0;
  lim3d.y.max = 0;

  lim3d.z.enable = false; // Limit the atom positions in the interval in the z-direction.
  lim3d.z.min = 0;
  lim3d.z.max = 0;

  // =========================================================
  // Check if command-line arguments are provided
  // =========================================================

  if (argc > 1)
  {
    // Parse command-line arguments
    for (int i = 1; i < argc; ++i)
    {
      std::string arg = argv[i];
      if (arg == "-f" && i + 1 < argc)
      {
        xyz_filename = argv[++i];
      }
      else if (arg == "-s" && i + 3 < argc)
      {
        start.x = std::atof(argv[++i]);
        start.y = std::atof(argv[++i]);
        start.z = std::atof(argv[++i]);
      }
      else if (arg == "-d" && i + 3 < argc)
      {
        direction.x = std::atof(argv[++i]);
        direction.y = std::atof(argv[++i]);
        direction.z = std::atof(argv[++i]);
      }
      else if (arg == "-l" && i + 1 < argc)
      {
        totalLength = std::atof(argv[++i]);
      }
      else if (arg == "-b" && i + 1 < argc)
      {
        numBins = std::atoi(argv[++i]);
      }
      else if (arg == "-start" && i + 1 < argc)
      {
        atomTypeStart = std::atoi(argv[++i]);
      }
      else if (arg == "-end" && i + 1 < argc)
      {
        atomTypeEnd = std::atoi(argv[++i]);
      }
      else if (arg == "-limx" && i + 1 < argc)
      {
        lim3d.x.enable = true;
        lim3d.x.min = std::atof(argv[++i]);
        lim3d.x.max = std::atof(argv[++i]);
      }
      else if (arg == "-limy" && i + 1 < argc)
      {
        lim3d.y.enable = true;
        lim3d.y.min = std::atof(argv[++i]);
        lim3d.y.max = std::atof(argv[++i]);
      }
      else if (arg == "-limz" && i + 1 < argc)
      {
        lim3d.z.enable = true;
        lim3d.z.min = std::atof(argv[++i]);
        lim3d.z.max = std::atof(argv[++i]);
      }
    }
  }

  // =========================================================
  // Generating Output file name
  // =========================================================

  std::string xyz_filename_without_xyz = "";
  std::string output_filename_prefix = "t_";
  {
    // Find the position of ".xyz" in the string
    size_t pos = xyz_filename.find(".xyz");

    // If ".xyz" exists in the string, remove it
    if (pos != std::string::npos)
    {
      // Remove ".xyz" from the string
      xyz_filename_without_xyz = xyz_filename.substr(0, pos);
    }
  }

  // Generate output folder name
  std::string outputFolder = generateOutputFolderName(xyz_filename_without_xyz);

  // Create output directory
  if (!directoryExists(outputFolder))
  {
    int result = mkdir(outputFolder.c_str(), 0777);
    if (result != 0)
    {
      std::cerr << "Error: Failed to create output directory." << std::endl;
      return 1;
    }
  }

  // =========================================================
  // Print running parameters
  // =========================================================

  std::string runningParameters = "Running Parameters:\n";
  runningParameters += "XYZ Filename: " + xyz_filename + "\n";
  runningParameters += "Start Position: (" + std::to_string(start.x) + ", " + std::to_string(start.y) + ", " + std::to_string(start.z) + ")\n";
  runningParameters += "Direction: (" + std::to_string(direction.x) + ", " + std::to_string(direction.y) + ", " + std::to_string(direction.z) + ")\n";
  runningParameters += "Total Length: " + std::to_string(totalLength) + "\n";
  runningParameters += "Number of Bins: " + std::to_string(numBins) + "\n";
  runningParameters += "Atom Type Start: " + std::to_string(atomTypeStart) + "\n";
  runningParameters += "Atom Type End: " + std::to_string(atomTypeEnd) + "\n";
  if (lim3d.x.enable)
    runningParameters += "Limited in x direction: [" + std::to_string(lim3d.x.min) + " , " + std::to_string(lim3d.x.max) + "]\n";
  if (lim3d.y.enable)
    runningParameters += "Limited in y direction: [" + std::to_string(lim3d.y.min) + " , " + std::to_string(lim3d.y.max) + "]\n";
  if (lim3d.z.enable)
    runningParameters += "Limited in z direction: [" + std::to_string(lim3d.z.min) + " , " + std::to_string(lim3d.z.max) + "]\n";

  // Print the running parameters
  std::cout << runningParameters << std::endl;

  // Write input parameters to file
  std::ofstream outputFile(outputFolder + "/input_parameters.txt");
  if (outputFile.is_open())
  {

    outputFile << runningParameters << std::endl;

    outputFile.close();
    std::cout << "Input parameters written to: " << outputFolder << "/input_parameters.txt" << std::endl;
  }
  else
  {
    std::cerr << "Error: Unable to open file for writing input parameters." << std::endl;
    return 1;
  }

  // =========================================================
  // Pre-processing
  // =========================================================

  double binLength = totalLength / numBins;
  if (direction.Normalize() != 0)
  {
    std::cout << "Program stopped due to an error \n";
    return 1;
  }

  // =========================================================
  // Importing XYZ file
  // =========================================================

  std::vector<Snapshot> snapshots = readXYZFile(xyz_filename);
  std::cout << "Total number of imported snapshots: " << snapshots.size() << "\n";

  // =========================================================
  // Calculating the linear distribution function and exporting it
  // =========================================================

  for (int atomType = atomTypeStart; atomType <= atomTypeEnd; atomType++)
  {
    std::cout << "atomType " << atomType << std::flush;

    std::vector<double> averageDistribution = calculateAverageDistribution(snapshots, atomType, start, direction, numBins, totalLength);

    std::string output_filename = output_filename_prefix + std::to_string(atomType);

    std::cout << " filename '" << output_filename << "'" << std::endl;

    std::ofstream ofs((outputFolder + "/" + output_filename).c_str());

    // Output average distribution
    for (int i = 0; i < numBins; ++i)
    {
      ofs << i << " " << binLength * i << " " << averageDistribution[i] << std::endl;
    }
  }

  std::cout << "Calculation finished. \n";

  return 0;
}
