
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

#ifndef CAVIAR_OBJECTS_FORCEFIELD_LJ_H
#define CAVIAR_OBJECTS_FORCEFIELD_LJ_H


#include "caviar/objects/force_field.h"
#include "caviar/utility/python_utils_dec.h"

namespace caviar {
namespace objects {
class Atom_data;
namespace atom_data {
class Basic;
}
namespace force_field {

/**
 * This class calculates LJ potential for the particles.
 * 
 */
class Lj : public Force_field {
public:
  Lj (class CAVIAR *);
  ~Lj () {};
  
  bool read (class caviar::interpreter::Parser *);
  void verify_settings ();
  void calculate_acceleration ();
public:
  bool input_by_array;
  std::vector<std::vector<Real_t>> epsilon,sigma;
  bool make_off_diagonal_vectors;




  bool input_by_atom;
  // epsilon - sigma of a single type. inter-type values will be deduced using these
  std::vector<Real_t> epsilon_atom, sigma_atom; 

//  std::vector< T > py_list_to_std_vector;
 // boost::python::list std_vector_to_py_list

  // some helper variables, used for debugging
  bool jump_fix, monitor_jump;
  double jump_tol;

  bool wca; //Week-Chandler-Anderson (WCA) potential activated.
  bool cutoff_list_activated;
  std::vector<std::vector<Real_t>> cutoff_list; // list of cutoffs when it is needed.
                                                // for example in WCA potentials

  FC_PYDEC_SETGET_PTR(atom_data,Atom_data);
  FC_PYDEC_SETGET_PTR(domain,Domain);
  FC_PYDEC_SETGET_PTR(neighborlist,Neighborlist);

  FC_PYDEC_SETGET_STDVEC2D(epsilon,Real_t);  
  FC_PYDEC_SETGET_STDVEC2D(sigma,Real_t);
  FC_PYDEC_SETGET_STDVEC(epsilon_atom,Real_t);  
  FC_PYDEC_SETGET_STDVEC(sigma_atom,Real_t);
  FC_PYDEC_SETGET_STDVEC2D(cutoff_list,Real_t);

};

void export_py_Lj ();

} //force_field
} //objects
} // namespace caviar

#endif
