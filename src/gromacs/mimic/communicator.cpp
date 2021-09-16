/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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

#include "communicator.h"

#include "config.h"

#include <unordered_set>

#include "gromacs/math/units.h"
#include "gromacs/utility/arrayref.h"

// Include headers from MiMiC library
#include <DataTypes.h>
#include <MessageApi.h>

namespace gmx
{

void MimicCommunicator::init()
{
    char path[GMX_PATH_MAX];
    gmx_getcwd(path, GMX_PATH_MAX);
    MCL_init_client(path);
}

void MimicCommunicator::sendInitData(gmx_mtop_t* mtop, ArrayRef<const RVec> coords)
{
    MCL_send(&mtop->natoms, 1, TYPE_INT, 0);
    MCL_send(&mtop->atomtypes.nr, 1, TYPE_INT, 0);

    std::vector<int>        atomTypes;
    std::vector<int>        nAtomsMol;
    std::vector<int>        idOrder;
    std::vector<double>     charges;
    std::vector<double>     masses(mtop->atomtypes.nr, -1);
    std::vector<int>        elements(mtop->atomtypes.nr, -1);
    std::vector<int>        bonds;
    std::vector<double>     bondLengths;
    std::unordered_set<int> existingTypes;

    atomTypes.reserve(static_cast<size_t>(mtop->natoms));
    charges.reserve(static_cast<size_t>(mtop->natoms));

    int offset = 0;
    for (const gmx_molblock_t& molblock : mtop->molblock)
    {
        gmx_moltype_t* type = &mtop->moltype[molblock.type];
        for (int mol = 0; mol < molblock.nmol; ++mol)
        {
            int nconstr  = type->ilist[F_CONSTR].size() / 3;
            int nconstrc = type->ilist[F_CONSTRNC].size() / 3;
            int nsettle  = type->ilist[F_SETTLE].size() / 4;

            for (int ncon = 0; ncon < nconstr + nconstrc; ++ncon)
            {
                int contype = type->ilist[F_CONSTR].iatoms[0];
                int at1     = type->ilist[F_CONSTR].iatoms[1];
                int at2     = type->ilist[F_CONSTR].iatoms[2];
                bonds.push_back(offset + at1 + 1);
                bonds.push_back(offset + at2 + 1);
                bondLengths.push_back(static_cast<double>(mtop->ffparams.iparams[contype].constr.dA)
                                      / c_bohr2Nm);
            }

            for (int ncon = 0; ncon < nsettle; ++ncon)
            {
                t_iatom ox;
                t_iatom h1;
                t_iatom h2;

                int contype = type->ilist[F_SETTLE].iatoms[0];

                ox = type->ilist[F_SETTLE].iatoms[1];
                h1 = type->ilist[F_SETTLE].iatoms[2];
                h2 = type->ilist[F_SETTLE].iatoms[3];

                bonds.push_back(offset + ox + 1);
                bonds.push_back(offset + h1 + 1);

                bonds.push_back(offset + ox + 1);
                bonds.push_back(offset + h2 + 1);

                bonds.push_back(offset + h1 + 1);
                bonds.push_back(offset + h2 + 1);
                bondLengths.push_back(static_cast<double>(mtop->ffparams.iparams[contype].constr.dA)
                                      / gmx::c_bohr2Nm);
                bondLengths.push_back(static_cast<double>(mtop->ffparams.iparams[contype].constr.dA)
                                      / gmx::c_bohr2Nm);
                bondLengths.push_back(static_cast<double>(mtop->ffparams.iparams[contype].constr.dB)
                                      / gmx::c_bohr2Nm);
            }

            nAtomsMol.push_back(type->atoms.nr);
            for (int at = 0; at < type->atoms.nr; ++at)
            {
                int  atomtype = type->atoms.atom[at].type;
                auto charge   = static_cast<double>(type->atoms.atom[at].q);
                idOrder.push_back(offset + 1);
                offset++;
                atomTypes.push_back(atomtype + 1);
                charges.push_back(charge);
                if (existingTypes.insert(atomtype).second)
                {
                    masses[atomtype]   = type->atoms.atom[at].m;
                    elements[atomtype] = type->atoms.atom[at].atomnumber;
                }
            }
        }
    }
    // sending atom types
    MCL_send(&*atomTypes.begin(), mtop->natoms, TYPE_INT, 0);

    int max_multipole_order = 0;
    // sending multipole orders
    MCL_send(&max_multipole_order, 1, TYPE_INT, 0);

    int nMolecules = nAtomsMol.size();
    // sending molecule number
    MCL_send(&nMolecules, 1, TYPE_INT, 0);

    // sending number of atoms per molecules
    MCL_send(&*nAtomsMol.begin(), nAtomsMol.size(), TYPE_INT, 0);

    int nBonds = bonds.size() / 2;
    // sending number of bond constraints
    MCL_send(&nBonds, 1, TYPE_INT, 0);

    // sending number of angle constraints
    MCL_send(&max_multipole_order, 1, TYPE_INT, 0);

    if (nBonds > 0)
    {
        // sending bonded atoms indices
        MCL_send(&*bonds.begin(), bonds.size(), TYPE_INT, 0);

        // sending bond lengths
        MCL_send(&*bondLengths.begin(), bondLengths.size(), TYPE_DOUBLE, 0);
    }

    // sending array of atomic charges
    MCL_send(&*charges.begin(), mtop->natoms, TYPE_DOUBLE, 0);

    // sending array of atomic masses
    MCL_send(&*masses.begin(), mtop->atomtypes.nr, TYPE_DOUBLE, 0);

    // sending ids of atoms per molecule
    MCL_send(&*idOrder.begin(), idOrder.size(), TYPE_INT, 0);

    // sending list of elements
    MCL_send(&*elements.begin(), mtop->atomtypes.nr, TYPE_INT, 0);

    std::vector<double> convertedCoords;
    for (auto& coord : coords)
    {
        convertedCoords.push_back(static_cast<double>(coord[0]) / gmx::c_bohr2Nm);
        convertedCoords.push_back(static_cast<double>(coord[1]) / gmx::c_bohr2Nm);
        convertedCoords.push_back(static_cast<double>(coord[2]) / gmx::c_bohr2Nm);
    }

    // sending array of coordinates
    MCL_send(&*convertedCoords.begin(), 3 * mtop->natoms, TYPE_DOUBLE, 0);
}

int64_t MimicCommunicator::getStepNumber()
{
    int steps;
    MCL_receive(&steps, 1, TYPE_INT, 0);
    return steps;
}

void MimicCommunicator::getCoords(ArrayRef<RVec> x, const int natoms)
{
    std::vector<double> coords(natoms * 3);
    MCL_receive(&*coords.begin(), 3 * natoms, TYPE_DOUBLE, 0);
    for (int j = 0; j < natoms; ++j)
    {
        x[j][0] = static_cast<real>(coords[j * 3] * gmx::c_bohr2Nm);
        x[j][1] = static_cast<real>(coords[j * 3 + 1] * gmx::c_bohr2Nm);
        x[j][2] = static_cast<real>(coords[j * 3 + 2] * gmx::c_bohr2Nm);
    }
}

void MimicCommunicator::sendEnergies(real energy)
{
    double convertedEnergy = energy / (gmx::c_hartree2Kj * gmx::c_avogadro);
    MCL_send(&convertedEnergy, 1, TYPE_DOUBLE, 0);
}

void MimicCommunicator::sendForces(ArrayRef<RVec> forces, int natoms)
{
    std::vector<double> convertedForce;
    for (int j = 0; j < natoms; ++j)
    {
        convertedForce.push_back(static_cast<real>(forces[j][0]) / gmx::c_hartreeBohr2Md);
        convertedForce.push_back(static_cast<real>(forces[j][1]) / gmx::c_hartreeBohr2Md);
        convertedForce.push_back(static_cast<real>(forces[j][2]) / gmx::c_hartreeBohr2Md);
    }
    MCL_send(&*convertedForce.begin(), convertedForce.size(), TYPE_DOUBLE, 0);
}

void MimicCommunicator::finalize()
{
    MCL_destroy();
}

} // namespace gmx
