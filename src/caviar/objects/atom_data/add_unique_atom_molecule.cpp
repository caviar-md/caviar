
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
        for (unsigned int j = 0; j < molecule_struct_owned[mi_1].atomic_bond_vector.size(); ++j)
        {
            auto b = molecule_struct_owned[mi_1].atomic_bond_vector[j];

            if (FC_COMPARE_PAIRS(b.id_1, b.id_2, bond.id_1, bond.id_2))
            {
                molecule_struct_owned[mi_1].atomic_bond_vector.erase(molecule_struct_owned[mi_1].atomic_bond_vector.begin() + j);
                remove_atom_id_from_molecule(b.id_1, mi_1);
                remove_atom_id_from_molecule(b.id_2, mi_1);
                bond_found = true;
                break;
            }
        }
        if (!bond_found)
        {
            error->all(FC_FILE_LINE_FUNC, "Can't remove bond. It does not exist");
            return;
        }

        for (unsigned int j = 0; j < molecule_struct_owned[mi_1].atomic_angle_vector.size(); ++j)
        {
            auto a = molecule_struct_owned[mi_1].atomic_angle_vector[j];
            if (FC_COMPARE_PAIRS(bond.id_1, bond.id_2, a.id_1, a.id_2) ||
                FC_COMPARE_PAIRS(bond.id_1, bond.id_2, a.id_2, a.id_3))
            {
                molecule_struct_owned[mi_1].atomic_angle_vector.erase(molecule_struct_owned[mi_1].atomic_angle_vector.begin() + j);
                remove_atom_id_from_molecule(a.id_1, mi_1);
                remove_atom_id_from_molecule(a.id_2, mi_1);
                remove_atom_id_from_molecule(a.id_3, mi_1);
                // break; no break! It can be shared between two angles
            }
        }

        for (unsigned int j = 0; j < molecule_struct_owned[mi_1].atomic_properdihedral_vector.size(); ++j)
        {
            auto p = molecule_struct_owned[mi_1].atomic_properdihedral_vector[j];
            if (FC_COMPARE_PAIRS(bond.id_1, bond.id_2, p.id_1, p.id_2) ||
                FC_COMPARE_PAIRS(bond.id_1, bond.id_2, p.id_2, p.id_3) ||
                FC_COMPARE_PAIRS(bond.id_1, bond.id_2, p.id_3, p.id_4))
            {
                molecule_struct_owned[mi_1].atomic_properdihedral_vector.erase(molecule_struct_owned[mi_1].atomic_properdihedral_vector.begin() + j);
                remove_atom_id_from_molecule(p.id_1, mi_1);
                remove_atom_id_from_molecule(p.id_2, mi_1);
                remove_atom_id_from_molecule(p.id_3, mi_1);
                remove_atom_id_from_molecule(p.id_4, mi_1);
                // break; no break! It can be shared between three proper dihedrals
            }
        }

        //------------------------------------------
        // remove mono-atomic  molecules  //
        //------------------------------------------
        if (molecule_struct_owned[mi_1].atomic_bond_vector.size() == 0)
        {
            num_molecules--;
            molecule_struct_owned.erase(molecule_struct_owned.begin() + mi_1);
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
        for (unsigned int j = 0; j < molecule_struct_owned[mi_1].atomic_angle_vector.size(); ++j)
        {
            auto a = molecule_struct_owned[mi_1].atomic_angle_vector[j];
            if (
                (FC_COMPARE_PAIRS(angle.id_1, angle.id_2, a.id_1, a.id_2) &&
                 FC_COMPARE_PAIRS(angle.id_2, angle.id_3, a.id_2, a.id_3)) ||
                (FC_COMPARE_PAIRS(angle.id_1, angle.id_2, a.id_2, a.id_3) &&
                 FC_COMPARE_PAIRS(angle.id_2, angle.id_3, a.id_1, a.id_2)))
            {
                molecule_struct_owned[mi_1].atomic_angle_vector.erase(molecule_struct_owned[mi_1].atomic_angle_vector.begin() + j);
                remove_atom_id_from_molecule(angle.id_1, mi_1);
                remove_atom_id_from_molecule(angle.id_2, mi_1);
                remove_atom_id_from_molecule(angle.id_3, mi_1);
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
        for (unsigned int j = 0; j < molecule_struct_owned[mi_1].atomic_properdihedral_vector.size(); ++j)
        {
            auto p = molecule_struct_owned[mi_1].atomic_properdihedral_vector[j];
            if (
                (FC_COMPARE_PAIRS(pd.id_1, pd.id_2, p.id_1, p.id_2) &&
                 FC_COMPARE_PAIRS(pd.id_3, pd.id_4, p.id_3, p.id_4)) ||
                (FC_COMPARE_PAIRS(pd.id_1, pd.id_2, p.id_3, p.id_4) &&
                 FC_COMPARE_PAIRS(pd.id_3, pd.id_4, p.id_1, p.id_2)))
            {
                molecule_struct_owned[mi_1].atomic_properdihedral_vector.erase(molecule_struct_owned[mi_1].atomic_properdihedral_vector.begin() + j);
                remove_atom_id_from_molecule(pd.id_1, mi_1);
                remove_atom_id_from_molecule(pd.id_2, mi_1);
                remove_atom_id_from_molecule(pd.id_3, mi_1);
                remove_atom_id_from_molecule(pd.id_4, mi_1);
                
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
    int bindex_1 = atom_id_to_index[bond.id_1];
    int bindex_2 = atom_id_to_index[bond.id_2];
    int molecule_index_1 = atom_struct_owned.molecule_index[bindex_1];
    int molecule_index_2 = atom_struct_owned.molecule_index[bindex_2];

    atom_struct_owned.atomic_bond_count[bindex_1]++;
    atom_struct_owned.atomic_bond_count[bindex_2]++;

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
        molecule_struct_owned.resize(num_molecules);

        int new_molecule_index = num_molecules - 1;

        //------------------------------------------
        // add the bond
        //------------------------------------------
        molecule_struct_owned[new_molecule_index].atomic_bond_vector.emplace_back(bond);
        atom_struct_owned.molecule_index[bindex_1] = new_molecule_index;
        atom_struct_owned.molecule_index[bindex_2] = new_molecule_index;
    }
    else if (molecule_index_1 == -1)
    {
        molecule_struct_owned[molecule_index_2].atomic_bond_vector.push_back(bond);
        atom_struct_owned.molecule_index[bindex_1] = molecule_index_1;
    }
    else if (molecule_index_2 == -1)
    {
        molecule_struct_owned[molecule_index_1].atomic_bond_vector.push_back(bond);
        atom_struct_owned.molecule_index[bindex_2] = molecule_index_1;
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

        auto vs = molecule_struct_owned[molecule_index_higher].atomic_bond_vector.size();

        if (vs > 0)
        {
            molecule_struct_owned[molecule_index_lower].atomic_bond_vector.reserve(vs);

            for (const auto &i : molecule_struct_owned[molecule_index_higher].atomic_bond_vector)
            {
                molecule_struct_owned[molecule_index_lower].atomic_bond_vector.emplace_back(i);
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_1]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_2]] = molecule_index_lower;
                add_atom_id_to_molecule(i.id_1, molecule_index_lower);
                add_atom_id_to_molecule(i.id_2, molecule_index_lower);
                remove_atom_id_from_molecule(i.id_1, molecule_index_higher);
                remove_atom_id_from_molecule(i.id_2, molecule_index_higher);
            }
            molecule_struct_owned[molecule_index_higher].atomic_bond_vector.clear();

        }

        molecule_struct_owned[molecule_index_lower].atomic_bond_vector.emplace_back(bond);

        //---------------------------------------------------------
        // After two molecules are merged by a new bond, existing
        // angles and properdihedrals have to be merged too.
        //--------------------------------------------------------
        vs = molecule_struct_owned[molecule_index_higher].atomic_angle_vector.size();
        if (vs > 0)
        {
            molecule_struct_owned[molecule_index_lower].atomic_angle_vector.reserve(vs);

            for (const auto &i : molecule_struct_owned[molecule_index_higher].atomic_angle_vector)
            {
                molecule_struct_owned[molecule_index_lower].atomic_angle_vector.emplace_back(i);
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_1]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_2]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_3]] = molecule_index_lower;

                add_atom_id_to_molecule(i.id_1, molecule_index_lower);
                add_atom_id_to_molecule(i.id_2, molecule_index_lower);
                add_atom_id_to_molecule(i.id_3, molecule_index_lower);
                remove_atom_id_from_molecule(i.id_1, molecule_index_higher);
                remove_atom_id_from_molecule(i.id_2, molecule_index_higher);           
                remove_atom_id_from_molecule(i.id_3, molecule_index_higher);
            }
            molecule_struct_owned[molecule_index_higher].atomic_angle_vector.clear();

        }

        vs = molecule_struct_owned[molecule_index_higher].atomic_properdihedral_vector.size();
        if (vs > 0)
        {
            molecule_struct_owned[molecule_index_lower].atomic_properdihedral_vector.reserve(vs);

            for (const auto &i : molecule_struct_owned[molecule_index_higher].atomic_properdihedral_vector)
            {
                molecule_struct_owned[molecule_index_lower].atomic_properdihedral_vector.emplace_back(i);
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_1]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_2]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_3]] = molecule_index_lower;
                atom_struct_owned.molecule_index[atom_id_to_index[i.id_4]] = molecule_index_lower;

                add_atom_id_to_molecule(i.id_1, molecule_index_lower);
                add_atom_id_to_molecule(i.id_2, molecule_index_lower);
                add_atom_id_to_molecule(i.id_3, molecule_index_lower);
                add_atom_id_to_molecule(i.id_4, molecule_index_lower);
                remove_atom_id_from_molecule(i.id_1, molecule_index_higher);
                remove_atom_id_from_molecule(i.id_2, molecule_index_higher);
                remove_atom_id_from_molecule(i.id_3, molecule_index_higher);
                remove_atom_id_from_molecule(i.id_4, molecule_index_higher);
            }

            molecule_struct_owned[molecule_index_higher].atomic_properdihedral_vector.clear();

        }
    }
    else if ((molecule_index_1 == molecule_index_2) && (molecule_index_1 != -1))
    {
        molecule_struct_owned[molecule_index_1].atomic_bond_vector.push_back(bond);
    }
}

void Atom_data::add_atomic_angle(const atom_data::Angle &a)
{

    int mi_1 = atom_struct_owned.molecule_index[atom_id_to_index[a.id_1]];
    int mi_2 = atom_struct_owned.molecule_index[atom_id_to_index[a.id_2]];
    int mi_3 = atom_struct_owned.molecule_index[atom_id_to_index[a.id_3]];

    if ((mi_1 == mi_2) && (mi_2 == mi_3))
    {
        molecule_struct_owned[mi_1].atomic_angle_vector.emplace_back(a);
        add_atom_id_to_molecule(a.id_1, mi_1);
        add_atom_id_to_molecule(a.id_2, mi_1);
        add_atom_id_to_molecule(a.id_3, mi_1);
    }
    else
    {
        error->all(FC_FILE_LINE_FUNC, "New atomic angles must be in the same molecule");
    }
}

void Atom_data::add_atomic_properdihedral(const atom_data::Proper_dihedral &pd)
{

    int mi_1 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_1]];
    int mi_2 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_2]];
    int mi_3 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_3]];
    int mi_4 = atom_struct_owned.molecule_index[atom_id_to_index[pd.id_4]];

    if ((mi_1 == mi_2) && (mi_2 == mi_3) && (mi_3 == mi_4))
    {
        molecule_struct_owned[mi_1].atomic_properdihedral_vector.emplace_back(pd);      
        add_atom_id_to_molecule(pd.id_1, mi_1);
        add_atom_id_to_molecule(pd.id_2, mi_1);
        add_atom_id_to_molecule(pd.id_3, mi_1);
        add_atom_id_to_molecule(pd.id_4, mi_1);
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
        for (unsigned int j = 0; j < molecule_struct_owned[mi_1].atomic_bond_vector.size(); ++j)
        {
            auto b = molecule_struct_owned[mi_1].atomic_bond_vector[j];
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
        for (unsigned int j = 0; j < molecule_struct_owned[mi_1].atomic_angle_vector.size(); ++j)
        {
            auto a = molecule_struct_owned[mi_1].atomic_angle_vector[j];
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
        for (unsigned int j = 0; j < molecule_struct_owned[mi_1].atomic_properdihedral_vector.size(); ++j)
        {
            auto p = molecule_struct_owned[mi_1].atomic_properdihedral_vector[j];
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
    molecule_struct_owned.resize(num_molecules);
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

        add_atom_id_to_molecule(id, new_molecule_index);
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

        molecule_struct_owned[new_molecule_index].atomic_bond_vector.emplace_back(dummy_atomic_bond);

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

        molecule_struct_owned[new_molecule_index].atomic_angle_vector.emplace_back(dummy_atomic_angle);
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

        molecule_struct_owned[new_molecule_index].atomic_properdihedral_vector.emplace_back(dummy_atomic_properdihedral);
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



void Atom_data::add_atom_id_to_molecule(int atom_id, int molecule_index)
{
    std::vector<int>::iterator position = std::find(molecule_struct_owned[molecule_index].atom_list.begin(), molecule_struct_owned[molecule_index].atom_list.end(), atom_id);
    if ( position != molecule_struct_owned[molecule_index].atom_list.end() )        
    {
        molecule_struct_owned[molecule_index].atom_list.emplace_back(atom_id);
    }
}


void Atom_data::remove_atom_id_from_molecule(int atom_id, int molecule_index)
{
    std::vector<int>::iterator position = std::find(molecule_struct_owned[molecule_index].atom_list.begin(), molecule_struct_owned[molecule_index].atom_list.end(), atom_id);
    if (position != molecule_struct_owned[molecule_index].atom_list.end()) // == myVector.end() means the element was not found
        molecule_struct_owned[molecule_index].atom_list.erase(position);
    //else
    //    error->all(FC_FILE_LINE_FUNC, "Can't remove atomic from molecule because it does not exist in the molecule::Atom_list.");

}

CAVIAR_NAMESPACE_CLOSE
