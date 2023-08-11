
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
#include "caviar/objects/unique/atom.h"
#include "caviar/objects/unique/atom_group.h"
#include "caviar/objects/unique/atom_list.h"
#include "caviar/objects/unique/molecule.h"
#include "caviar/objects/unique/molecule_group.h"
#include "caviar/objects/unique/molecule_list.h"

#include <algorithm>
#include <random>

CAVIAR_NAMESPACE_OPEN

bool Atom_data::add_atom(caviar::unique::Atom &a)
{

    const auto id = get_num_of_atoms_global();
    const auto t = a.type;
    const auto pos = a.pos_tot();
    const auto vel = a.vel_tot();
    if (position_inside_local_domain(pos))
        return add_atom(id, t, pos, vel);
    return false;
}

bool Atom_data::add_atom(caviar::unique::Atom_group &ag)
{
    for (auto &&a : ag.atoms)
        add_atom(a);
    return true;
}

bool Atom_data::add_atom(caviar::unique::Atom_list &al)
{
    for (auto &&a : al.atoms)
        add_atom(*a);
    return true;
}

inline bool FC_COMPARE_PAIRS(int a1, int a2, int b1, int b2)
{
    if ((a1 == b1 && a2 == b2) || (a1 == b2 && a2 == b1))
        return true;
    return false;
}

void Atom_data::remove_atomic_bond(const atom_data::Bond &bond)
{

    int mi_1 = atom_struct_owned.molecule_index[atom_id_to_index[bond.id_1]];
    int mi_2 = atom_struct_owned.molecule_index[atom_id_to_index[bond.id_2]];
    if (mi_1 == -1 || mi_2 == -1)
        return;

    if (mi_1 == mi_2)
    {
        bool bond_found = false;
        for (unsigned int j = 0; j < molecule_struct_owned.atomic_bond_vector[mi_1].size(); ++j)
        {
            auto b = molecule_struct_owned.atomic_bond_vector[mi_1][j];

            if (FC_COMPARE_PAIRS(b.id_1, b.id_2, bond.id_1, bond.id_2))
            {
                molecule_struct_owned.atomic_bond_vector[mi_1].erase(molecule_struct_owned.atomic_bond_vector[mi_1].begin() + j);

                bond_found = true;
                break;
            }
        }
        if (!bond_found)
        {
            error->all(FC_FILE_LINE_FUNC, "Can't remove bond. It does not exist");
            return;
        }

        for (unsigned int j = 0; j < molecule_struct_owned.atomic_angle_vector[mi_1].size(); ++j)
        {
            auto a = molecule_struct_owned.atomic_angle_vector[mi_1][j];
            if (FC_COMPARE_PAIRS(bond.id_1, bond.id_2, a.id_1, a.id_2) ||
                FC_COMPARE_PAIRS(bond.id_1, bond.id_2, a.id_2, a.id_3))
            {
                molecule_struct_owned.atomic_angle_vector[mi_1].erase(molecule_struct_owned.atomic_angle_vector[mi_1].begin() + j);
                // break; no break! It can be shared between two angles
            }
        }

        for (unsigned int j = 0; j < molecule_struct_owned.atomic_properdihedral_vector[mi_1].size(); ++j)
        {
            auto p = molecule_struct_owned.atomic_properdihedral_vector[mi_1][j];
            if (FC_COMPARE_PAIRS(bond.id_1, bond.id_2, p.id_1, p.id_2) ||
                FC_COMPARE_PAIRS(bond.id_1, bond.id_2, p.id_2, p.id_3) ||
                FC_COMPARE_PAIRS(bond.id_1, bond.id_2, p.id_3, p.id_4))
            {
                molecule_struct_owned.atomic_properdihedral_vector[mi_1].erase(molecule_struct_owned.atomic_properdihedral_vector[mi_1].begin() + j);
                // break; no break! It can be shared between three proper dihedrals
            }
        }

        //------------------------------------------
        // remove mono-atomic  molecules  //
        //------------------------------------------
        if (molecule_struct_owned.atomic_bond_vector[mi_1].size() == 0)
        {
            num_molecules--;
            molecule_struct_owned.atomic_bond_vector.erase(molecule_struct_owned.atomic_bond_vector.begin() + mi_1);
            molecule_struct_owned.atomic_angle_vector.erase(molecule_struct_owned.atomic_angle_vector.begin() + mi_1);
            molecule_struct_owned.atomic_properdihedral_vector.erase(molecule_struct_owned.atomic_properdihedral_vector.begin() + mi_1);
        }
    }
    else
    {
        error->all(FC_FILE_LINE_FUNC, "Can't remove atomic bond because the described atoms are of different molecules.");
    }

    atom_struct_owned.atomic_bond_count[atom_id_to_index[bond.id_1]]--;
    atom_struct_owned.atomic_bond_count[atom_id_to_index[bond.id_2]]--;
}

void Atom_data::remove_atomic_angle(const atom_data::Angle &angle)
{

    int mi_1 = atom_struct_owned.molecule_index[atom_id_to_index[angle.id_1]];
    int mi_2 = atom_struct_owned.molecule_index[atom_id_to_index[angle.id_2]];
    int mi_3 = atom_struct_owned.molecule_index[atom_id_to_index[angle.id_3]];

    if ((mi_1 == mi_2) && (mi_2 == mi_3))
    {
        for (unsigned int j = 0; j < molecule_struct_owned.atomic_angle_vector[mi_1].size(); ++j)
        {
            auto a = molecule_struct_owned.atomic_angle_vector[mi_1][j];
            if (
                (FC_COMPARE_PAIRS(angle.id_1, angle.id_2, a.id_1, a.id_2) &&
                 FC_COMPARE_PAIRS(angle.id_2, angle.id_3, a.id_2, a.id_3)) ||
                (FC_COMPARE_PAIRS(angle.id_1, angle.id_2, a.id_2, a.id_3) &&
                 FC_COMPARE_PAIRS(angle.id_2, angle.id_3, a.id_1, a.id_2)))
            {
                molecule_struct_owned.atomic_angle_vector[mi_1].erase(molecule_struct_owned.atomic_angle_vector[mi_1].begin() + j);
                break;
            }
        }
    }
    else
    {
        error->all(FC_FILE_LINE_FUNC, "Can't remove atomic angle because the described atoms are of different molecules.");
    }
}

void Atom_data::remove_atomic_properdihedral(const atom_data::Proper_dihedral &pd)
{

    int mi_1 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_1]];
    int mi_2 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_2]];
    int mi_3 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_3]];
    int mi_4 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_4]];

    if ((mi_1 == mi_2) && (mi_2 == mi_3) && (mi_3 == mi_4))
    {
        for (unsigned int j = 0; j < molecule_struct_owned.atomic_properdihedral_vector[mi_1].size(); ++j)
        {
            auto p = molecule_struct_owned.atomic_properdihedral_vector[mi_1][j];
            if (
                (FC_COMPARE_PAIRS(pd.id_1, pd.id_2, p.id_1, p.id_2) &&
                 FC_COMPARE_PAIRS(pd.id_3, pd.id_4, p.id_3, p.id_4)) ||
                (FC_COMPARE_PAIRS(pd.id_1, pd.id_2, p.id_3, p.id_4) &&
                 FC_COMPARE_PAIRS(pd.id_3, pd.id_4, p.id_1, p.id_2)))
            {
                molecule_struct_owned.atomic_properdihedral_vector[mi_1].erase(molecule_struct_owned.atomic_properdihedral_vector[mi_1].begin() + j);
                break;
            }
        }
    }
    else
    {
        error->all(FC_FILE_LINE_FUNC, "Can't remove atomic proper dihedral because the described atoms are of different molecules.");
    }
}

void Atom_data::add_atomic_bond(const atom_data::Bond &bond)
{

    int molecule_index_1 = atom_struct_owned.molecule_index[atom_id_to_index[bond.id_1]];
    int molecule_index_2 = atom_struct_owned.molecule_index[atom_id_to_index[bond.id_2]];

    atom_struct_owned.atomic_bond_count[atom_id_to_index[bond.id_1]]++;
    atom_struct_owned.atomic_bond_count[atom_id_to_index[bond.id_2]]++;

    //---------------------------------------------------------
    // deducing whether a new molecule must be created or
    // choosing the molecule with the lower index to be the one
    // that the other atoms/molecule are merged into.
    //---------------------------------------------------------

    if (molecule_index_1 == -1 && molecule_index_2 == -1)
    {

        //------------------------------------------
        // increasing Molecules' containers size  //
        //------------------------------------------
        num_molecules++;
        molecule_struct_owned.atomic_bond_vector.resize(num_molecules);
        molecule_struct_owned.atomic_angle_vector.resize(num_molecules);
        molecule_struct_owned.atomic_properdihedral_vector.resize(num_molecules);
        int new_molecule_index = num_molecules - 1;

        //------------------------------------------
        // add the bond
        //------------------------------------------
        molecule_struct_owned.atomic_bond_vector[new_molecule_index].emplace_back(bond);
        atom_struct_owned.molecule_index[atom_id_to_index[bond.id_1]] = new_molecule_index;
        atom_struct_owned.molecule_index[atom_id_to_index[bond.id_2]] = new_molecule_index;
    }
    else if (molecule_index_1 == -1)
    {
        molecule_struct_owned.atomic_bond_vector[molecule_index_2].push_back(bond);
        atom_struct_owned.molecule_index[atom_id_to_index[bond.id_1]] = molecule_index_1;
    }
    else if (molecule_index_2 == -1)
    {
        molecule_struct_owned.atomic_bond_vector[molecule_index_1].push_back(bond);
        atom_struct_owned.molecule_index[atom_id_to_index[bond.id_2]] = molecule_index_1;
    }
    else if (molecule_index_1 != molecule_index_2)
    {
        int molecule_index_lower, molecule_index_higher;
        if (molecule_index_1 < molecule_index_2)
        {
            molecule_index_lower = molecule_index_1;
            molecule_index_higher = molecule_index_2;
        }
        else
        {
            molecule_index_lower = molecule_index_2;
            molecule_index_higher = molecule_index_1;
        }

        auto vs = molecule_struct_owned.atomic_bond_vector[molecule_index_higher].size();

        if (vs > 0)
        {
            molecule_struct_owned.atomic_bond_vector[molecule_index_lower].reserve(vs);

            for (const auto &i : molecule_struct_owned.atomic_bond_vector[molecule_index_higher])
            {
                molecule_struct_owned.atomic_bond_vector[molecule_index_lower].emplace_back(i);
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_1]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_2]] = molecule_index_lower;
            }
            molecule_struct_owned.atomic_bond_vector[molecule_index_higher].clear();
        }

        molecule_struct_owned.atomic_bond_vector[molecule_index_lower].emplace_back(bond);

        //---------------------------------------------------------
        // After two molecules are merged by a new bond, existing
        // angles and properdihedrals have to be merged too.
        //--------------------------------------------------------
        vs = molecule_struct_owned.atomic_angle_vector[molecule_index_higher].size();
        if (vs > 0)
        {
            molecule_struct_owned.atomic_angle_vector[molecule_index_lower].reserve(vs);

            for (const auto &i : molecule_struct_owned.atomic_angle_vector[molecule_index_higher])
            {
                molecule_struct_owned.atomic_angle_vector[molecule_index_lower].emplace_back(i);
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_1]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_2]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_3]] = molecule_index_lower;
            }
            molecule_struct_owned.atomic_angle_vector[molecule_index_higher].clear();
        }

        vs = molecule_struct_owned.atomic_properdihedral_vector[molecule_index_higher].size();
        if (vs > 0)
        {
            molecule_struct_owned.atomic_properdihedral_vector[molecule_index_lower].reserve(vs);

            for (const auto &i : molecule_struct_owned.atomic_properdihedral_vector[molecule_index_higher])
            {
                molecule_struct_owned.atomic_properdihedral_vector[molecule_index_lower].emplace_back(i);
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_1]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_2]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_3]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_4]] = molecule_index_lower;
            }

            molecule_struct_owned.atomic_properdihedral_vector[molecule_index_higher].clear();
        }
    }
    else if ((molecule_index_1 == molecule_index_2) && (molecule_index_1 != -1))
    {
        molecule_struct_owned.atomic_bond_vector[molecule_index_1].push_back(bond);
    }
}

void Atom_data::add_atomic_angle(const atom_data::Angle &angle)
{

    int mi_1 = atom_struct_owned.molecule_index[atom_id_to_index[angle.id_1]];
    int mi_2 = atom_struct_owned.molecule_index[atom_id_to_index[angle.id_2]];
    int mi_3 = atom_struct_owned.molecule_index[atom_id_to_index[angle.id_3]];

    if ((mi_1 == mi_2) && (mi_2 == mi_3))
    {
        molecule_struct_owned.atomic_angle_vector[mi_1].emplace_back(angle);
    }
    else
    {
        error->all(FC_FILE_LINE_FUNC, "New atomic angles must be in the same molecule");
    }
}

void Atom_data::add_atomic_properdihedral(const atom_data::Proper_dihedral &proper_dihedral)
{

    int mi_1 = atom_struct_owned.molecule_index[atom_id_to_index[proper_dihedral.id_1]];
    int mi_2 = atom_struct_owned.molecule_index[atom_id_to_index[proper_dihedral.id_2]];
    int mi_3 = atom_struct_owned.molecule_index[atom_id_to_index[proper_dihedral.id_3]];
    int mi_4 = atom_struct_owned.molecule_index[atom_id_to_index[proper_dihedral.id_4]];

    if ((mi_1 == mi_2) && (mi_2 == mi_3) && (mi_3 == mi_4))
    {
        molecule_struct_owned.atomic_properdihedral_vector[mi_1].emplace_back(proper_dihedral);
    }
    else
    {
        error->all(FC_FILE_LINE_FUNC, "New atomic proper dihedrals must be in the same molecule");
    }
}

bool Atom_data::check_atomic_bond_exist(const atom_data::Bond &bond)
{
    int mi_1 = atom_struct_owned.molecule_index[atom_id_to_index[bond.id_1]];
    int mi_2 = atom_struct_owned.molecule_index[atom_id_to_index[bond.id_2]];
    if (mi_1 == -1 || mi_2 == -1)
        return false;
    if (mi_1 == mi_2)
    {
        for (unsigned int j = 0; j < molecule_struct_owned.atomic_bond_vector[mi_1].size(); ++j)
        {
            auto b = molecule_struct_owned.atomic_bond_vector[mi_1][j];
            if (FC_COMPARE_PAIRS(b.id_1, b.id_2, bond.id_1, bond.id_2))
            {
                return true;
            }
        }
    }
    return false;
}

bool Atom_data::check_atomic_angle_exist(const atom_data::Angle &angle)
{

    int mi_1 = atom_struct_owned.molecule_index[atom_id_to_index[angle.id_1]];
    int mi_2 = atom_struct_owned.molecule_index[atom_id_to_index[angle.id_2]];
    int mi_3 = atom_struct_owned.molecule_index[atom_id_to_index[angle.id_3]];

    if ((mi_1 == mi_2) && (mi_2 == mi_3))
    {
        for (unsigned int j = 0; j < molecule_struct_owned.atomic_angle_vector[mi_1].size(); ++j)
        {
            auto a = molecule_struct_owned.atomic_angle_vector[mi_1][j];
            if (
                (FC_COMPARE_PAIRS(angle.id_1, angle.id_2, a.id_1, a.id_2) &&
                 FC_COMPARE_PAIRS(angle.id_2, angle.id_3, a.id_2, a.id_3)) ||
                (FC_COMPARE_PAIRS(angle.id_1, angle.id_2, a.id_2, a.id_3) &&
                 FC_COMPARE_PAIRS(angle.id_2, angle.id_3, a.id_1, a.id_2)))
            {
                return true;
            }
        }
    }
    return false;
}

bool Atom_data::check_atomic_properdihedral_exist(const atom_data::Proper_dihedral &pd)
{

    int mi_1 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_1]];
    int mi_2 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_2]];
    int mi_3 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_3]];
    int mi_4 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_4]];

    if ((mi_1 == mi_2) && (mi_2 == mi_3) && (mi_3 == mi_4))
    {
        for (unsigned int j = 0; j < molecule_struct_owned.atomic_properdihedral_vector[mi_1].size(); ++j)
        {
            auto p = molecule_struct_owned.atomic_properdihedral_vector[mi_1][j];
            if (
                (FC_COMPARE_PAIRS(pd.id_1, pd.id_2, p.id_1, p.id_2) &&
                 FC_COMPARE_PAIRS(pd.id_3, pd.id_4, p.id_3, p.id_4)) ||
                (FC_COMPARE_PAIRS(pd.id_1, pd.id_2, p.id_3, p.id_4) &&
                 FC_COMPARE_PAIRS(pd.id_3, pd.id_4, p.id_1, p.id_2)))
            {
                return true;
            }
        }
    }
    return false;
}

bool Atom_data::add_molecule(caviar::unique::Molecule &m)
{

    //------------------------------------------
    // increasing Molecules' containers size  //
    //------------------------------------------
    num_molecules++;
    molecule_struct_owned.atomic_bond_vector.resize(num_molecules);
    molecule_struct_owned.atomic_angle_vector.resize(num_molecules);
    molecule_struct_owned.atomic_properdihedral_vector.resize(num_molecules);
    int new_molecule_index = num_molecules - 1;

    //-------------------------------------
    // extract the atoms of the Molecule //
    //-------------------------------------
    std::vector<int> types;
    std::vector<Vector<double>> pos, vel;

    m.extract_all_e_pos_vel(types, pos, vel);

    //--------------------------------------------------------------
    // check whether the center of mass of atoms of the molecule
    // is inside the MPI domain or not. If not, ignore the molecule
    //--------------------------------------------------------------
    Vector<double> pos_cm(0, 0, 0);
    for (unsigned int i = 0; i < pos.size(); ++i)
    {
        pos_cm += pos[i];
    }
    pos_cm *= 1.0 / pos.size();

    if (!position_inside_local_domain(pos_cm))
        return false;

    //----------------------------------------------------------
    // add a copy of molecule's atom into Atom data as new atoms
    // making a map of the local indices to global
    //----------------------------------------------------------

    std::vector<int> indices(pos.size(), -1);
    std::vector<int> ids(pos.size(), -1);

    for (unsigned int i = 0; i < pos.size(); ++i)
    {
        const auto id = get_num_of_atoms_global();
        const auto index = atom_struct_owned.position.size();
        add_atom(id, types[i], pos[i], vel[i]);
        indices[i] = index;
        ids[i] = id;
        atom_struct_owned.molecule_index[index] = new_molecule_index;
    }

    //---------------------------
    // Adding the atomic bonds //
    //---------------------------
    for (unsigned int j = 0; j < m.atomic_bond.size(); ++j)
    {
        auto dummy_atomic_bond = m.atomic_bond[j];
        // Changes the local indices to global ones
        dummy_atomic_bond.id_1 = ids[dummy_atomic_bond.id_1];//indices[dummy_atomic_bond.id_1];
        dummy_atomic_bond.id_2 = ids[dummy_atomic_bond.id_2];//indices[dummy_atomic_bond.id_2];

        molecule_struct_owned.atomic_bond_vector[new_molecule_index].emplace_back(dummy_atomic_bond);

        atom_struct_owned.atomic_bond_count[dummy_atomic_bond.id_1]++;
        atom_struct_owned.atomic_bond_count[dummy_atomic_bond.id_2]++;
    }

    //----------------------------
    // Adding the atomic angles //
    //----------------------------
    for (unsigned int j = 0; j < m.atomic_angle.size(); ++j)
    {
        auto dummy_atomic_angle = m.atomic_angle[j];
        // Changes the local indices to global ones
        dummy_atomic_angle.id_1 = ids[dummy_atomic_angle.id_1];
        dummy_atomic_angle.id_2 = ids[dummy_atomic_angle.id_2];
        dummy_atomic_angle.id_3 = ids[dummy_atomic_angle.id_3];

        molecule_struct_owned.atomic_angle_vector[new_molecule_index].emplace_back(dummy_atomic_angle);
    }

    //--------------------------------------
    // Adding the atomic proper dihedrals //
    //--------------------------------------
    for (unsigned int j = 0; j < m.atomic_properdihedral.size(); ++j)
    {
        auto dummy_atomic_properdihedral = m.atomic_properdihedral[j];
        // Changes the local indices to global ones
        dummy_atomic_properdihedral.id_1 = ids[dummy_atomic_properdihedral.id_1];
        dummy_atomic_properdihedral.id_2 = ids[dummy_atomic_properdihedral.id_2];
        dummy_atomic_properdihedral.id_3 = ids[dummy_atomic_properdihedral.id_3];
        dummy_atomic_properdihedral.id_4 = ids[dummy_atomic_properdihedral.id_4];

        molecule_struct_owned.atomic_properdihedral_vector[new_molecule_index].emplace_back(dummy_atomic_properdihedral);
    }

    return true;
}

bool Atom_data::add_molecule(caviar::unique::Molecule_group &mg)
{
    for (auto &&a : mg.molecules)
        add_molecule(a);
    return true;
}

bool Atom_data::add_molecule(caviar::unique::Molecule_list &ml)
{
    for (auto &&a : ml.molecules)
        add_molecule(*a);
    return true;
}

CAVIAR_NAMESPACE_CLOSE
