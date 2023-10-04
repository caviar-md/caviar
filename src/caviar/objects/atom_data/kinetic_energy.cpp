
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
#include "caviar/interpreter/output.h"
#include "caviar/objects/domain.h"
#include "caviar/utility/common_template_functions.h"
#include "caviar/objects/unique/time_function_3d.h"
#include <algorithm>
#ifdef CAVIAR_WITH_MPI
#include <mpi.h>
#endif
CAVIAR_NAMESPACE_OPEN

double Atom_data::kinetic_energy_mpi_domain(const int t)
{

  double e_owned = 0.0;
  if (velocity_offset == nullptr)
  {

    switch (get_n_r_df())
    {
    case 0:
    {

      for (unsigned int i = 0; i < atom_struct_owned.position.size(); ++i)
      {

        if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
          continue;

        if (t > -1 && t != static_cast<int>(atom_struct_owned.type[i]))
          continue; // KINETIC ENERGY OF A TYPE
        e_owned += atom_type_params.mass[atom_struct_owned.type[i]] * (atom_struct_owned.velocity[i] * atom_struct_owned.velocity[i]);
      }
    }
    break;

    case (1):
    case (2):
    case (3):
    {

      auto v_cm = owned_velocity_cm();
      for (unsigned int i = 0; i < atom_struct_owned.position.size(); ++i)
      {

        if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
          continue;

        if (t > -1 && t != static_cast<int>(atom_struct_owned.type[i]))
          continue; // KINETIC ENERGY OF A TYPE
        auto v_i = atom_struct_owned.velocity[i] - v_cm;
        e_owned += atom_type_params.mass[atom_struct_owned.type[i]] * (v_i * v_i);
      }
    }
    break;

    case (4):
    case (5):
    case (6):
    {

      auto v_cm = owned_velocity_cm_mpi_domain();
      auto p_cm = owned_position_cm_mpi_domain();

      auto L_cm = owned_angular_momentum_cm_mpi_domain();
      auto I_cm = owned_inertia_tensor_cm_mpi_domain();

      std::array<std::array<double, 3>, 3> I_cm_inverse = {{{0.0, 0.0, 0.0},
                                                            {0.0, 0.0, 0.0},
                                                            {0.0, 0.0, 0.0}}};
      Vector<Real_t> I_i_L(0, 0, 0);
      bool correct_result = false;
      if (matrix_inverse_3d(I_cm, I_cm_inverse) != 0)
      {
        matrix_Vector_product_3d(I_cm_inverse, L_cm, I_i_L);
        correct_result = true;
      }
      else
      {
        output->warning("Error in matrix_inverse_3d for Atom_data::kinetic_energy.");
        // error->all(FC_FILE_LINE_FUNC, "Error in Inertia Tensor Inverse Calculations for Atom_data::kinetic_energy.");
      }

      if (correct_result)
      {

        for (unsigned int i = 0; i < atom_struct_owned.position.size(); ++i)
        {

          if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
            continue;

          if (t > -1 && t != static_cast<int>(atom_struct_owned.type[i]))
            continue; // KINETIC ENERGY OF A TYPE

          auto v_i = atom_struct_owned.velocity[i] - v_cm - (cross_product(I_i_L, atom_struct_owned.position[i] - p_cm));
          e_owned += atom_type_params.mass[atom_struct_owned.type[i]] * (v_i * v_i);
        }
      }
    }
    break;
    }
  }
  else
  { // TODO: mix this block into ' switch (get_n_r_df())' cases
    for (unsigned int i = 0; i < atom_struct_owned.position.size(); ++i)
    {
      if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;

      if (t > -1 && t != static_cast<int>(atom_struct_owned.type[i]))
        continue; // KINETIC ENERGY OF A TYPE
      auto fixed_vel = atom_struct_owned.velocity[i] + velocity_offset->current_value;
      e_owned += atom_type_params.mass[atom_struct_owned.type[i]] * (fixed_vel * fixed_vel);
    }
  }
  e_owned *= 0.5;

  return e_owned;
}

double Atom_data::kinetic_energy(const int t)
{
  double e_local = kinetic_energy_mpi_domain(t);

  double e_total = 0.0;

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  e_total = e_local;
#elif defined(CAVIAR_WITH_MPI)
  MPI_Allreduce(&e_local, &e_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  e_total = e_local;
#endif
  return e_total;
}

void Atom_data::add_to_temperature(double, int)
{
  // pressure_mpi_domain_ += v;
  error->all(FC_FILE_LINE_FUNC, "Not implemented");
}

void Atom_data::finalize_temperature()
{

  if (!temperature_process)
    return;

  reset_temperature();

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  finalize_temperature_mpi_domain();
  temperature_ = temperature_mpi_domain_;
#elif defined(CAVIAR_WITH_MPI)
  finalize_temperature_total();
  finalize_temperature_mpi_domain();
#else
  finalize_temperature_mpi_domain();
  temperature_ = temperature_mpi_domain_;
#endif
}

void Atom_data::finalize_temperature_total()
{
  if (k_b < 0)
    error->all(FC_FILE_LINE_FUNC, "k_b is not set.");

  // by equipartition theorem, k = 1/2 * k_b * N_df * T.
  // It is implemented according to
  // 'Berendsen and Nose-Hoover thermostats Victor Ruhle August 8, 2007'
  temperature_ = 2.0 * kinetic_energy() / (k_b * degree_of_freedoms());
}

void Atom_data::finalize_temperature_mpi_domain()
{
  if (k_b < 0)
    error->all(FC_FILE_LINE_FUNC, "k_b is not set.");

  // by equipartition theorem, k = 1/2 * k_b * N_df * T.
  // It is implemented according to
  // 'Berendsen and Nose-Hoover thermostats Victor Ruhle August 8, 2007'
  temperature_mpi_domain_ = 2.0 * kinetic_energy_mpi_domain() / (k_b * degree_of_freedoms_mpi_domain());
}

void Atom_data::add_to_pressure(double, int)
{
  // pressure_mpi_domain_ += v;
  error->all(FC_FILE_LINE_FUNC, "Not implemented");
}

void Atom_data::finalize_pressure()
{

  if (!pressure_process)
    return;

  reset_pressure();

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  finalize_pressure_mpi_domain();
  pressure_ = pressure_mpi_domain_;
#elif defined(CAVIAR_WITH_MPI)
  finalize_pressure_total();
  finalize_pressure_mpi_domain();
#else
  finalize_pressure_mpi_domain();
  pressure_ = pressure_mpi_domain_;
#endif
}

void Atom_data::finalize_pressure_mpi_domain()
{
  auto &acc = atom_struct_owned.acceleration;

  auto d_diff = (domain->upper_local - domain->lower_local);
  double volume = d_diff.x * d_diff.y * d_diff.z;

  double p2 = 0;
  int64_t pos_size = atom_struct_owned.position.size();

  double p2_old = 0;
  double Num_active = 0;

  // XXX note that this STATIC value may affect NPT ensembles
  caviar::Vector<double> domain_dx = {(domain->upper_global.x - domain->lower_global.x),
                                      (domain->upper_global.y - domain->lower_global.y),
                                      (domain->upper_global.z - domain->lower_global.z)};

  bool finish = false;
  for (int i = 0; i < pos_size; ++i)
  {
#ifdef CAVIAR_WITH_MPI
    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;
#endif
    Num_active++;

    int type = atom_struct_owned.type[i];
    double mass = atom_type_params.mass[type];
    auto f2 = acc[i] * mass;


    double fix_x =atom_struct_owned.msd_domain_cross[i].x * domain_dx.x;
    double fix_y =atom_struct_owned.msd_domain_cross[i].y * domain_dx.y;
    double fix_z =atom_struct_owned.msd_domain_cross[i].z * domain_dx.z; 

    Vector<double> pos_msd_i = {atom_struct_owned.position[i].x + fix_x,
                                atom_struct_owned.position[i].y + fix_y,
                                atom_struct_owned.position[i].z + fix_z};
    if (fix_x >0 || fix_y>0 || fix_z >0)
    {
      std::cout << "i: " << i <<" fix: "<< fix_x <<" , " << fix_y << " , " << fix_z << std::endl;
      finish = true;
    }
    // p2_old += (f2) * atom_struct_owned.position[i];
    p2_old += (f2)*pos_msd_i;
  }
  p2_old = p2_old / (3.0 * volume);

  p2 = (virialConstraint + virialExternalForce + virialForce) / (3.0 * volume);
  //p2 = (0 + virialExternalForce + virialForce) / (3.0 * volume); // XXXXXXXX

  //std::cout << "Constraint: " << virialConstraint / (3.0 * volume) << " ExternalForce:" << virialExternalForce / (3.0 * volume) << " Force: " << virialForce / (3.0 * volume) << std::endl;

  double p1 = (Num_active * k_b * temperature_mpi_domain()) / volume;
    // std::cout << "x p2_old: " << p2_old << std::endl;

  pressure_mpi_domain_ = p1 + p2;
  //pressure_mpi_domain_ = p1 + p2_old;
  //std::cout << "p1: " << p1 << " p2:" << p2 << " p2_old: " << p2_old << std::endl;
  // std::cout << "pressure: " << pressure_mpi_domain_ << std::endl;
  if (pressure_mpi_domain_ < 0)
    output->warning("Negative pressure");
  if (finish)
    error->all(FC_FILE_LINE_FUNC, "finish");
}

void Atom_data::finalize_pressure_total()
{

  auto d_diff = (domain->upper_global - domain->lower_global);
  double volume = d_diff.x * d_diff.y * d_diff.z;

  double p2 = 0;
  int64_t pos_size = atom_struct_owned.position.size();

  double Num_active = 0;
  for (int i = 0; i < pos_size; ++i)
  {
#ifdef CAVIAR_WITH_MPI
    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;
#endif
    Num_active++;

    int type = atom_struct_owned.type[i];
    double mass = atom_type_params.mass[type];
    auto f = atom_struct_owned.acceleration[i] * mass;
    p2 += f * atom_struct_owned.position[i];
  }

  double p2_total = 0;
#ifdef CAVIAR_WITH_MPI

  MPI_Allreduce(&p2, &p2_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  p2_total = p2;
#endif

  double p1 = (Num_active * k_b * temperature()) / volume;

  p2_total = p2_total / (3.0 * volume);

  pressure_ = p1 + p2_total;
}

Vector<Real_t> Atom_data::owned_position_cm()
{
  double p_cm_f[3] = {0.0, 0.0, 0.0};
  double mass_sum = 0.0;
  double mass_sum_total = 0.0;
  auto p_size = atom_struct_owned.position.size(); // MPI check
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ : p_cm, mass_sum)
#endif
  for (unsigned int i = 0; i < p_size; ++i)
  {
    auto type_i = atom_struct_owned.type[i];
    auto mass_i = atom_type_params.mass[type_i];
    mass_sum_total += mass_i;

    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;

    mass_sum += mass_i;
    p_cm_f[0] += atom_struct_owned.position[i].x * mass_i;
    p_cm_f[1] += atom_struct_owned.position[i].y * mass_i;
    p_cm_f[2] += atom_struct_owned.position[i].z * mass_i;
  }

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#elif defined(CAVIAR_WITH_MPI)
  double p_cm_f_total[3] = {0.0, 0.0, 0.0};
  MPI_Allreduce(&p_cm_f, &p_cm_f_total, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  p_cm_f[0] = p_cm_f_total[0];
  p_cm_f[1] = p_cm_f_total[1];
  p_cm_f[2] = p_cm_f_total[2];

  mass_sum = mass_sum_total;
#else

#endif
  Vector<Real_t> p_cm{p_cm_f[0], p_cm_f[1], p_cm_f[2]};

  p_cm = p_cm / mass_sum;
  return p_cm;
}

Vector<Real_t> Atom_data::owned_position_cm_mpi_domain()
{
  Vector<Real_t> p_cm{0.0, 0.0, 0.0};
  double mass_sum = 0.0;
  auto p_size = atom_struct_owned.position.size(); // MPI check
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ : p_cm, mass_sum)
#endif
  for (unsigned int i = 0; i < p_size; ++i)
  {
    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;

    auto type_i = atom_struct_owned.type[i];
    auto mass_i = atom_type_params.mass[type_i];
    mass_sum += mass_i;
    p_cm += atom_struct_owned.position[i] * mass_i;
  }
  p_cm = p_cm / mass_sum;
  return p_cm;
}

Vector<Real_t> Atom_data::owned_velocity_cm()
{
  double v_cm_f[3] = {0.0, 0.0, 0.0};

  double mass_sum = 0.0;
  double mass_sum_total = 0.0;

  auto p_size = atom_struct_owned.velocity.size(); // MPI check
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ : v_cm, mass_sum)
#endif
  for (unsigned int i = 0; i < p_size; ++i)
  {
    auto type_i = atom_struct_owned.type[i];
    auto mass_i = atom_type_params.mass[type_i];

    mass_sum_total += mass_i;

    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;

    mass_sum += mass_i;
    v_cm_f[0] += atom_struct_owned.velocity[i].x * mass_i;
    v_cm_f[1] += atom_struct_owned.velocity[i].y * mass_i;
    v_cm_f[2] += atom_struct_owned.velocity[i].z * mass_i;
  }

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#elif defined(CAVIAR_WITH_MPI)
  double v_cm_f_total[3] = {0.0, 0.0, 0.0};
  MPI_Allreduce(&v_cm_f, &v_cm_f_total, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  v_cm_f[0] = v_cm_f_total[0];
  v_cm_f[1] = v_cm_f_total[1];
  v_cm_f[2] = v_cm_f_total[2];

  mass_sum = mass_sum_total;
#else

#endif
  Vector<Real_t> v_cm{v_cm_f[0], v_cm_f[1], v_cm_f[2]};

  v_cm = v_cm / mass_sum;
  return v_cm;
}

Vector<Real_t> Atom_data::owned_velocity_cm_mpi_domain()
{
  Vector<Real_t> v_cm{0.0, 0.0, 0.0};
  double mass_sum = 0.0;
  auto p_size = atom_struct_owned.velocity.size(); // MPI check
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ : v_cm, mass_sum)
#endif
  for (unsigned int i = 0; i < p_size; ++i)
  {
    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;

    auto type_i = atom_struct_owned.type[i];
    auto mass_i = atom_type_params.mass[type_i];
    mass_sum += mass_i;
    v_cm += atom_struct_owned.velocity[i] * mass_i;
  }
  v_cm = v_cm / mass_sum;
  return v_cm;
}

Vector<double> Atom_data::owned_angular_momentum_cm(const Vector<double> &p_cm)
{

  double L_cm_f[3] = {0.0, 0.0, 0.0};

  auto p_size = atom_struct_owned.position.size(); // MPI check
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ : L_cm)
#endif
  for (unsigned int i = 0; i < p_size; ++i)
  {
    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;

    auto type_i = atom_struct_owned.type[i];
    auto mass_i = atom_type_params.mass[type_i];
    auto res = mass_i * cross_product(atom_struct_owned.position[i] - p_cm, atom_struct_owned.velocity[i]);
    L_cm_f[0] += res.x;
    L_cm_f[1] += res.y;
    L_cm_f[2] += res.z;
  }

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#elif defined(CAVIAR_WITH_MPI)
  double L_cm_f_total[3] = {0.0, 0.0, 0.0};
  MPI_Allreduce(&L_cm_f, &L_cm_f_total, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  L_cm_f[0] = L_cm_f_total[0];
  L_cm_f[1] = L_cm_f_total[1];
  L_cm_f[2] = L_cm_f_total[2];

#else

#endif
  Vector<double> L_cm{L_cm_f[0], L_cm_f[1], L_cm_f[2]};

  return L_cm;
}

Vector<double> Atom_data::owned_angular_momentum_cm_mpi_domain(const Vector<double> &p_cm)
{

  Vector<double> L_cm{0.0, 0.0, 0.0};

  auto p_size = atom_struct_owned.position.size(); // MPI check
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ : L_cm)
#endif
  for (unsigned int i = 0; i < p_size; ++i)
  {
    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;

    auto type_i = atom_struct_owned.type[i];
    auto mass_i = atom_type_params.mass[type_i];

    L_cm += mass_i * cross_product(atom_struct_owned.position[i] - p_cm, atom_struct_owned.velocity[i]);
  }

  return L_cm;
}

std::array<std::array<double, 3>, 3> Atom_data::owned_inertia_tensor_cm(const Vector<double> &p_cm)
{

  double I_cm_flat[9] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0};

  auto p_size = atom_struct_owned.position.size(); // MPI check
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ : it_cm)
#endif
  for (unsigned int i = 0; i < p_size; ++i)
  {
    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;

    auto type_i = atom_struct_owned.type[i];
    auto mass_i = atom_type_params.mass[type_i];
    auto p = atom_struct_owned.position[i] - p_cm;    // relative position
    I_cm_flat[0] += mass_i * (p.y * p.y + p.z * p.z); //=I_cm[0][0]
    I_cm_flat[1] += mass_i * (p.x * p.x + p.z * p.z); //=I_cm[1][1]
    I_cm_flat[2] += mass_i * (p.x * p.x + p.y * p.y); //=I_cm[2][2]

    I_cm_flat[3] -= mass_i * (p.x * p.y); //=I_cm[0][1]
    I_cm_flat[4] -= mass_i * (p.y * p.z); //=I_cm[1][2]
    I_cm_flat[5] -= mass_i * (p.z * p.x); //=I_cm[2][0]
  }
  I_cm_flat[6] = I_cm_flat[3]; // I_cm[0][1]=I_cm[1][0]
  I_cm_flat[7] = I_cm_flat[4]; // I_cm[1][2]=I_cm[2][1]
  I_cm_flat[8] = I_cm_flat[5]; // I_cm[2][0]=I_cm[0][2]

  std::array<std::array<double, 3>, 3> I_cm;

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#elif defined(CAVIAR_WITH_MPI)
  double I_cm_flat_total[9] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0};
  MPI_Allreduce(&I_cm_flat, &I_cm_flat_total, 9, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  I_cm_flat[0] = I_cm_flat_total[0];
  I_cm_flat[1] = I_cm_flat_total[1];
  I_cm_flat[2] = I_cm_flat_total[2];
  I_cm_flat[3] = I_cm_flat_total[3];
  I_cm_flat[4] = I_cm_flat_total[4];
  I_cm_flat[5] = I_cm_flat_total[5];
  I_cm_flat[6] = I_cm_flat_total[6];
  I_cm_flat[7] = I_cm_flat_total[7];
  I_cm_flat[8] = I_cm_flat_total[8];

#else

#endif
  I_cm[0][0] = I_cm_flat[0];
  I_cm[1][1] = I_cm_flat[1];
  I_cm[2][2] = I_cm_flat[2];
  I_cm[0][1] = I_cm_flat[3];
  I_cm[1][2] = I_cm_flat[4];
  I_cm[2][0] = I_cm_flat[5];
  I_cm[1][0] = I_cm_flat[6];
  I_cm[2][1] = I_cm_flat[7];
  I_cm[0][2] = I_cm_flat[8];

  return I_cm;
}

std::array<std::array<double, 3>, 3> Atom_data::owned_inertia_tensor_cm_mpi_domain(const Vector<double> &p_cm)
{

  std::array<std::array<double, 3>, 3> I_cm = {{{0.0, 0.0, 0.0},
                                                {0.0, 0.0, 0.0},
                                                {0.0, 0.0, 0.0}}};

  auto p_size = atom_struct_owned.position.size(); // MPI check
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for reduction(+ : it_cm)
#endif
  for (unsigned int i = 0; i < p_size; ++i)
  {
    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;

    auto type_i = atom_struct_owned.type[i];
    auto mass_i = atom_type_params.mass[type_i];
    auto p = atom_struct_owned.position[i] - p_cm; // relative position
    I_cm[0][0] += mass_i * (p.y * p.y + p.z * p.z);
    I_cm[1][1] += mass_i * (p.x * p.x + p.z * p.z);
    I_cm[2][2] += mass_i * (p.x * p.x + p.y * p.y);

    I_cm[0][1] -= mass_i * (p.x * p.y);
    I_cm[1][2] -= mass_i * (p.y * p.z);
    I_cm[2][0] -= mass_i * (p.z * p.x);
  }

  I_cm[1][0] = I_cm[0][1];
  I_cm[2][1] = I_cm[1][2];
  I_cm[0][2] = I_cm[2][0];

  return I_cm;
}

int Atom_data::degree_of_freedoms()
{

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#elif defined(CAVIAR_WITH_MPI)

  int64_t dof_total = 0;
  int64_t dof_local = degree_of_freedoms_mpi_domain();
  MPI_Allreduce(&dof_local, &dof_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  return dof_total;

#else

#endif

  auto df = 3 * atom_struct_owned.position.size();

  auto sum_of_bonds = 0;
  auto sum_of_angles = 0;
  auto sum_of_dihedrals = 0;

  for (auto &&m : molecule_struct_owned)
  {
    sum_of_bonds += m.atomic_bond_vector.size();

    sum_of_angles += m.atomic_angle_vector.size();

    sum_of_dihedrals += m.atomic_properdihedral_vector.size();
  }

  return df - sum_of_bonds - sum_of_angles - sum_of_dihedrals - get_n_r_df();
}

int Atom_data::degree_of_freedoms_mpi_domain()
{

  int pos_size = atom_struct_owned.position.size();

  size_t df = 0;
  size_t sum_of_bonds = 0;
  size_t sum_of_angles = 0;
  size_t sum_of_dihedrals = 0;

  std::vector<uint8_t> molecule_selected(molecule_struct_owned.size() + 1, 0);

  for (int i = 0; i < pos_size; ++i)
  {
#ifdef CAVIAR_WITH_MPI
    if (atom_struct_owned.mpi_rank[i] != my_mpi_rank)
      continue;
#endif

    df += 3;

    auto mi = atom_struct_owned.molecule_index[i];

    if (mi == -1)
      continue;

    if (molecule_selected[mi] == 1)
      continue;

    molecule_selected[mi] = 1;

    sum_of_bonds += molecule_struct_owned[mi].atomic_bond_vector.size();

    sum_of_angles += molecule_struct_owned[mi].atomic_angle_vector.size();

    sum_of_dihedrals += molecule_struct_owned[mi].atomic_properdihedral_vector.size();
  }

  return df - sum_of_bonds - sum_of_angles - sum_of_dihedrals - get_n_r_df();
}

CAVIAR_NAMESPACE_CLOSE
