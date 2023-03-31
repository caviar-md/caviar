
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATICSPHERICALBOUNDARY_H
#define CAVIAR_OBJECTS_FORCEFIELD_ELECTROSTATICSPHERICALBOUNDARY_H

#include "caviar/objects/force_field.h"

CAVIAR_NAMESPACE_OPEN

namespace force_field {

/**
 * This class have electrostatic for charged particles in the spherical metalic
 * boundary. It has been used as a test for PLT method
 */
class Electrostatic_spherical_boundary : public Force_field {
public:
  Electrostatic_spherical_boundary (class CAVIAR *);
  ~Electrostatic_spherical_boundary () {};
  double potential (const Vector<double> &);
  double potential (const int);

  Vector<double> field (const Vector<double> &);
  Vector<double> field (const int);

  double energy();

  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();

  void calculate_image_charges();
  void initialize ();
  
  double potential_of_charge (const Vector<double> &, unsigned int i);
  double potential_of_image (const Vector<double> &, unsigned int i);
  
  double potential_of_charges (const Vector<double> &);// tot. potential at a point  due to charge only
  double potential_of_charges (unsigned int i);          // tot. potential at a point  due to charge only
  
  Vector<double> field_of_charge (const Vector<double> &, unsigned int i);
  Vector<double> field_of_image (const Vector<double> &, unsigned int i);
   
  Vector<double> field_of_charges (const Vector<double> &); // tot. field at a point due to real charges
  Vector<double> field_of_charges (unsigned int i); // tot. field on a charge from others

public:
  bool calculated_once;
  double k_electrostatic;
  Vector<double> external_field;

  double radius;         // radius of the spherical boundary
  Vector<double> center; // center of the spherical boundary
  double voltage;        // voltage on the spherical boundary 

  // this will optimise in the case of the existence of uncharged particles 
  // but it will cost one-to-one particle-image correspondance.
  bool uncharged_particles_optimization;

  struct {
    std::vector<Real_t> charge; // note that the charge can be different for
                                // each of images. So the charge in this struct
                                // is different from Atom_data's owned charge. 
    std::vector<Vector<Real_t>> position;
  } image;
 
};

} //force_field

CAVIAR_NAMESPACE_CLOSE

#endif
