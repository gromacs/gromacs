/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "MimicCommunicator.h"

#include <iostream>
#include <limits>
#include <vector>

#include <unordered_set>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/utility/fatalerror.h"

#ifdef GMX_MIMIC
#include <DataTypes.h>
#include <MessageApi.h>
#endif

void gmx::MimicCommunicator::init()
{
#ifdef GMX_MIMIC
    char path[GMX_PATH_MAX];
    gmx_getcwd(path, GMX_PATH_MAX);
    return MCL_init_client(path);
#else
    gmx_fatal(FARGS, "GROMACS is compiled without MiMiC support! Please, recompile with -DGMX_MIMIC=ON");
#endif
}

void gmx::MimicCommunicator::sendInitData(gmx_mtop_t gmx_unused *mtop, gmx::HostVector<gmx::RVec> gmx_unused coords)
{
#ifdef GMX_MIMIC
    MCL_send(&mtop->natoms, 1, TYPE_INT, 0);
    MCL_send(&mtop->atomtypes.nr, 1, TYPE_INT, 0);

    std::vector<int>        atomtypes;
    std::vector<int>        natoms_mol;
    std::vector<int>        id_order;
    std::vector<double>     charges;
    std::vector<double>     masses(mtop->atomtypes.nr, -1);
    std::vector<int>        elements(mtop->atomtypes.nr, -1);
    std::vector<int>        bonds;
    std::vector<double>     bond_lengths;
    std::unordered_set<int> exsiting_types;

    atomtypes.reserve((size_t)mtop->natoms);
    charges.reserve((size_t)mtop->natoms);

    int offset = 0;
    for (int i = 0; i < mtop->nmolblock; ++i)
    {
        gmx_molblock_t molblock = mtop->molblock[i];
        gmx_moltype_t  type     = mtop->moltype[molblock.type];
        for (int mol = 0; mol < molblock.nmol; ++mol)
        {
            int      nconstr  = type.ilist[F_CONSTR].nr / 3;
            int      nconstrc = type.ilist[F_CONSTRNC].nr / 3;
            int      nsettle  = type.ilist[F_SETTLE].nr / 4;

            t_iatom *at1 = mtop->moltype[mtop->molblock[i].type].ilist[F_CONSTR].iatoms;
            t_iatom *at2 = mtop->moltype[mtop->molblock[i].type].ilist[F_CONSTRNC].iatoms;

            for (int ncon = 0; ncon < nconstr + nconstrc; ++ncon)
            {
                t_iatom *atoms;
                atoms = constr_iatomptr(nconstr, at1, at2, ncon);
                int      contype = atoms[0];
                at1 = &atoms[1];
                at2 = &atoms[2];
                bonds.push_back((int &&)(offset + *at1 + 1));
                bonds.push_back((int &&)(offset + *at2 + 1));
                bond_lengths.push_back(static_cast<double>(mtop->ffparams.iparams[contype].constr.dA) / bohr);
            }

            for (int ncon = 0; ncon < nsettle; ++ncon)
            {
                t_iatom ox;
                t_iatom h1;
                t_iatom h2;

                int     contype = type.ilist[F_SETTLE].iatoms[0];

                ox = type.ilist[F_SETTLE].iatoms[1];
                h1 = type.ilist[F_SETTLE].iatoms[2];
                h2 = type.ilist[F_SETTLE].iatoms[3];

                bonds.push_back(offset + ox + 1);
                bonds.push_back(offset + h1 + 1);

                bonds.push_back(offset + ox + 1);
                bonds.push_back(offset + h2 + 1);

                bonds.push_back(offset + h1 + 1);
                bonds.push_back(offset + h2 + 1);
                bond_lengths.push_back(static_cast<double>(mtop->ffparams.iparams[contype].constr.dA) / bohr);
                bond_lengths.push_back(static_cast<double>(mtop->ffparams.iparams[contype].constr.dA) / bohr);
                bond_lengths.push_back(static_cast<double>(mtop->ffparams.iparams[contype].constr.dB) / bohr);
            }

            natoms_mol.push_back(type.atoms.nr);
            for (int at = 0; at < type.atoms.nr; ++at)
            {
                int  atomtype = type.atoms.atom[at].type;
                auto charge   = (double) type.atoms.atom[at].q;
                id_order.push_back(offset + 1);
                offset++;
                atomtypes.push_back(atomtype + 1);
                charges.push_back(charge);
                if (exsiting_types.insert(atomtype).second)
                {
                    masses[atomtype]   = type.atoms.atom[at].m;
                    elements[atomtype] = type.atoms.atom[at].atomnumber;
                }
            }
        }
    }
    // sending atom types
    MCL_send(&*atomtypes.begin(), mtop->natoms, TYPE_INT, 0);

    int max_multipole_order = 0;
    //sending multipole orders
    MCL_send(&max_multipole_order, 1, TYPE_INT, 0);

    int molecule_n = natoms_mol.size();
    // sending molecule number
    MCL_send(&molecule_n, 1, TYPE_INT, 0);

    // sending number of atoms per molecules
    MCL_send(&*natoms_mol.begin(), natoms_mol.size(), TYPE_INT, 0);

    int bonds_number = bonds.size() / 2;
    // sending number of bond constraints
    MCL_send(&bonds_number, 1, TYPE_INT, 0);

    // sending number of angle constraints
    MCL_send(&max_multipole_order, 1, TYPE_INT, 0);

    if (bonds_number > 0)
    {
        // sending bonded atoms indices
        MCL_send(&*bonds.begin(), bonds.size(), TYPE_INT, 0);

        // sending bond lengths
        MCL_send(&*bond_lengths.begin(), bond_lengths.size(), TYPE_DOUBLE, 0);
    }

    //sending array of atomic charges
    MCL_send(&*charges.begin(), mtop->natoms, TYPE_DOUBLE, 0);

    // sending array of atomic masses
    MCL_send(&*masses.begin(), mtop->atomtypes.nr, TYPE_DOUBLE, 0);

    // sending ids of atoms per molecule
    MCL_send(&*id_order.begin(), id_order.size(), TYPE_INT, 0);

    //sending list of elements
    MCL_send(&*elements.begin(), mtop->atomtypes.nr, TYPE_INT, 0);

    std::vector<double> trans_coord;
    for (int j = 0; j < coords.size(); ++j)
    {
        trans_coord.push_back(static_cast<double>(coords[j][0]) / bohr);
        trans_coord.push_back(static_cast<double>(coords[j][1]) / bohr);
        trans_coord.push_back(static_cast<double>(coords[j][2]) / bohr);
    }

    // sending array of coordinates
    MCL_send(&*trans_coord.begin(), 3 * mtop->natoms, TYPE_DOUBLE, 0);
#else
    gmx_fatal(FARGS, "GROMACS is compiled without MiMiC support! Please, recompile with -DGMX_MIMIC=ON");
#endif
}

void gmx::MimicCommunicator::getStepNumber(gmx_int64_t gmx_unused &nsteps)
{
#ifdef GMX_MIMIC
    int temp_step;
    MCL_receive(&temp_step, 1, TYPE_INT, 0);
    nsteps = temp_step;
#else
    gmx_fatal(FARGS, "GROMACS is compiled without MiMiC support! Please, recompile with -DGMX_MIMIC=ON");
#endif
}

void gmx::MimicCommunicator::getCoords(rvec gmx_unused *x, int gmx_unused natoms)
{
#ifdef GMX_MIMIC
    std::vector<double> trans_coords(natoms * 3);
    MCL_receive(&*trans_coords.begin(), 3 * natoms, TYPE_DOUBLE, 0);
    real                min_x = std::numeric_limits<real>::infinity();
    real                min_y = std::numeric_limits<real>::infinity();
    real                min_z = std::numeric_limits<real>::infinity();
    for (int j = 0; j < natoms; ++j)
    {
        x[j][0] = static_cast<real>(trans_coords[j * 3] * bohr);
        x[j][1] = static_cast<real>(trans_coords[j * 3 + 1] * bohr);
        x[j][2] = static_cast<real>(trans_coords[j * 3 + 2] * bohr);
        min_x   = std::min(x[j][0], min_x);
        min_y   = std::min(x[j][1], min_y);
        min_z   = std::min(x[j][2], min_z);
    }

    // We need to shift the origin of coordinates since CPMD sets {0,0,0} to be lower left corner of the box
    for (int j = 0; j < natoms; ++j)
    {
        x[j][0] -= min_x;
        x[j][1] -= min_y;
        x[j][2] -= min_z;
    }
#else
    gmx_fatal(FARGS, "GROMACS is compiled without MiMiC support! Please, recompile with -DGMX_MIMIC=ON");
#endif
}

void gmx::MimicCommunicator::sendEnergies(real gmx_unused energy)
{
#ifdef GMX_MIMIC
    double conv_energy = energy * hartree;
    MCL_send(&conv_energy, 1, TYPE_DOUBLE, 0);
#else
    gmx_fatal(FARGS, "GROMACS is compiled without MiMiC support! Please, recompile with -DGMX_MIMIC=ON");
#endif
}

void gmx::MimicCommunicator::sendForces(gmx::ArrayRef<gmx::RVec> gmx_unused forces, int gmx_unused natoms)
{
#ifdef GMX_MIMIC
    std::vector<double> trans_force;
    for (int j = 0; j < natoms; ++j)
    {
        trans_force.push_back(static_cast<real>(forces[j][0]) * hartree * bohr);
        trans_force.push_back(static_cast<real>(forces[j][1]) * hartree * bohr);
        trans_force.push_back(static_cast<real>(forces[j][2]) * hartree * bohr);
    }
    MCL_send(&*trans_force.begin(), trans_force.size(), TYPE_DOUBLE, 0);
#else
    gmx_fatal(FARGS, "GROMACS is compiled without MiMiC support! Please, recompile with -DGMX_MIMIC=ON");
#endif
}

void gmx::MimicCommunicator::finalize()
{
#ifdef GMX_MIMIC
    MCL_destroy();
#else
    gmx_fatal(FARGS, "GROMACS is compiled without MiMiC support! Please, recompile with -DGMX_MIMIC=ON");
#endif
}
