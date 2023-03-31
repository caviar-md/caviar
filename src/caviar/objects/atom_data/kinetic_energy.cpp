
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

#include "caviar/objects/atom_data.h"
#include "caviar/interpreter/communicator.h"
#include "caviar/interpreter/error.h"
#include "caviar/objects/domain.h"
#include "caviar/utility/common_template_functions.h"
#include "caviar/objects/unique/time_function_3d.h"
#include <algorithm>


namespace caviar {



double Atom_data::kinetic_energy () {
  double e_total = 0.0;

  double e_owned = 0.0;

  if (velocity_offset == nullptr)
  {
    switch (get_n_r_df()) 
    {
      case 0 : 
      {
        
        for (unsigned int i = 0; i < owned.position.size(); ++i) 
        {
          e_owned += owned.mass[owned.type[i]] * (owned.velocity[i] * owned.velocity[i]);
        }

      } break;

      case (1) : case (2): case (3): 
      {
        auto v_cm = owned_velocity_cm();
        for (unsigned int i = 0; i < owned.position.size(); ++i) 
        {
          auto v_i = owned.velocity[i] - v_cm;
          e_owned += owned.mass[owned.type[i]] * (v_i * v_i);
        }
      } break;

      case (4): case (5): case (6): 
      {
        
        auto v_cm = owned_velocity_cm();
        auto p_cm = owned_position_cm();          
        auto L_cm = owned_angular_momentum_cm();
        auto I_cm = owned_inertia_tensor_cm();
        
        std::vector<std::vector<Real_t> > I_cm_inverse (3, std::vector<Real_t> (3, 0.0)); 
        Vector<Real_t> I_i_L (0,0,0);
        bool correct_result = false;
        if (matrix_inverse(I_cm, I_cm_inverse)== 0) 
        {
            if (matrix_Vector_product(I_cm_inverse, L_cm, I_i_L) == 0)        
            {
                
                correct_result = true;
            }
        }
        
        if (correct_result)
        {
            for (unsigned int i = 0; i < owned.position.size(); ++i) 
            {
            
            auto v_i = owned.velocity[i] - v_cm - (cross_product(I_i_L, owned.position[i] - p_cm));
            e_owned += owned.mass[owned.type[i]] * (v_i * v_i);
            }
        }
        else
        {
            error->all(FC_FILE_LINE_FUNC,"Error in Inertia Tensor Inverse Calculations for Kinetic Energy...");    
        }
        
      } break;
    }
  }
  else
  { // TODO: mix this block into ' switch (get_n_r_df())' cases
      for (unsigned int i = 0; i < owned.position.size(); ++i) 
      {
        auto fixed_vel = owned.velocity[i] + velocity_offset->current_value;
        e_owned += owned.mass[owned.type[i]] * (fixed_vel * fixed_vel);
      }
  }
  e_owned *= 0.5;
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  e_total = e_owned;
#elif defined(CAVIAR_WITH_MPI)
  MPI_Allreduce (&e_owned, &e_total, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
#else
  e_total = e_owned;
#endif
  return e_total;
}


// KINETIC ENERGY OF A TYPE
double Atom_data::kinetic_energy (const int t) {
  
      
      
  double e_total = 0.0;

  double e_owned = 0.0;

  if (velocity_offset == nullptr)
  {
    switch (get_n_r_df()) 
    {
      case 0 : {
        for (unsigned int i = 0; i < owned.position.size(); ++i) 
        {
            if (t!=static_cast<int>(owned.type[i]))   continue; // KINETIC ENERGY OF A TYPE
            e_owned += owned.mass[owned.type[i]] * (owned.velocity[i] * owned.velocity[i]);
        }

      } break;

      case (1) : case (2): case (3): 
      {
        auto v_cm = owned_velocity_cm();
        for (unsigned int i = 0; i < owned.position.size(); ++i) 
        {
          if (t!=static_cast<int>(owned.type[i]))   continue; // KINETIC ENERGY OF A TYPE
          auto v_i = owned.velocity[i] - v_cm;
          e_owned += owned.mass[owned.type[i]] * (v_i * v_i);
        }
      } break;

      case (4): case (5): case (6):
      {
        
        auto v_cm = owned_velocity_cm();
        auto p_cm = owned_position_cm();    
        
        auto L_cm = owned_angular_momentum_cm();
        auto I_cm = owned_inertia_tensor_cm();
        
        std::vector<std::vector<Real_t> > I_cm_inverse (3, std::vector<Real_t> (3, 0.0)); 
        Vector<Real_t> I_i_L (0,0,0);
        bool correct_result = false;
        if (matrix_inverse(I_cm, I_cm_inverse)!= 0) 
        {
            if (matrix_Vector_product(I_cm_inverse, L_cm, I_i_L) != 0)
            {
                correct_result = true;
            }
        }
        
        if (correct_result)
        {
            for (unsigned int i = 0; i < owned.position.size(); ++i) 
            {
                if (t!=static_cast<int>(owned.type[i]))   continue; // KINETIC ENERGY OF A TYPE
            
                auto v_i = owned.velocity[i] - v_cm - (cross_product(I_i_L, owned.position[i] - p_cm));
                e_owned += owned.mass[owned.type[i]] * (v_i * v_i);
            }
        }
        else
        {
            error->all(FC_FILE_LINE_FUNC,"Error in Inertia Tensor Inverse Calculations for Kinetic Energy...");    
        }
        
      } break;
    }
  }
  else
  { // TODO: mix this block into ' switch (get_n_r_df())' cases
      for (unsigned int i = 0; i < owned.position.size(); ++i) 
      {
        if (t!=static_cast<int>(owned.type[i]))   continue; // KINETIC ENERGY OF A TYPE
        auto fixed_vel = owned.velocity[i] + velocity_offset->current_value;
        e_owned += owned.mass[owned.type[i]] * (fixed_vel * fixed_vel);
      }
  }

  e_owned *= 0.5;
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  e_total = e_owned;
#elif defined(CAVIAR_WITH_MPI)
  MPI_Allreduce (&e_owned, &e_total, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
#else
  e_total = e_owned;
#endif
  return e_total;

}

double Atom_data::temperature () {

  if (k_b < 0) error->all(FC_FILE_LINE_FUNC,"k_b is not set.");

  // by equipartition theorem, k = 1/2 * k_b * N_df * T.
  // It is implemented according to
  // 'Berendsen and Nose-Hoover thermostats Victor Ruhle August 8, 2007'
  return 2.0 * kinetic_energy() / (k_b * degree_of_freedoms());

}



} // namespace caviar


