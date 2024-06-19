/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#include "massrepartitioning.h"

#include <cstdio>

#include <algorithm>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "gromacs/fileio/warninp.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

static real smallestAtomMass(const gmx_mtop_t& mtop)
{
    real smallestMass = std::numeric_limits<real>::max();

    for (const gmx_moltype_t& moltype : mtop.moltype)
    {
        const t_atoms& atoms = moltype.atoms;

        for (int a = 0; a < atoms.nr; a++)
        {
            if (atoms.atom[a].m > 0)
            {
                smallestMass = std::min(smallestMass, atoms.atom[a].m);
            }
        }
    }

    return smallestMass;
}

void repartitionAtomMasses(gmx_mtop_t* mtop, const bool useFep, const real massFactor, WarningHandler* wi)
{
    const real smallestMass = smallestAtomMass(*mtop);

    const real minMass = massFactor * smallestMass;

    printf("The smallest mass in the system is %g, setting the minimum mass to %g\n", smallestMass, minMass);

    bool haveTooSmallPartnerMass    = false;
    int  numUnboundLightAtoms       = 0;
    int  numMultipleBoundLightAtoms = 0;

    for (gmx_moltype_t& moltype : mtop->moltype)
    {
        const t_atoms& atoms = moltype.atoms;

        std::vector<int> bondPartner(atoms.nr, -1);

        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            if ((interaction_function[ftype].flags & IF_CHEMBOND) == 0 || ftype == F_SETTLE)
            {
                continue;
            }

            const int c_numAtomsInBond = 2;

            GMX_RELEASE_ASSERT(NRAL(ftype) == c_numAtomsInBond,
                               "Expect only chemical bonds between two atoms");

            gmx::ArrayRef<const int> iatoms = moltype.ilist[ftype].iatoms;

            for (int i = 0; i < iatoms.ssize(); i += 1 + c_numAtomsInBond)
            {
                int& partner1 = bondPartner[iatoms[i + 1]];
                if (partner1 == -1)
                {
                    partner1 = iatoms[i + 2];
                }
                else
                {
                    partner1 = -2;
                }

                int& partner2 = bondPartner[iatoms[i + 2]];
                if (partner2 == -1)
                {
                    partner2 = iatoms[i + 1];
                }
                else
                {
                    partner2 = -2;
                }
            }
        }

        ArrayRef<const int> settleIatoms = moltype.ilist[F_SETTLE].iatoms;
        for (int i = 0; i < settleIatoms.ssize(); i += 1 + NRAL(F_SETTLE))
        {
            bondPartner[settleIatoms[i + 1]] = 0;

            int& partner2 = bondPartner[settleIatoms[i + 2]];
            if (partner2 == -1)
            {
                partner2 = settleIatoms[i + 1];
            }
            else
            {
                partner2 = -2;
            }

            int& partner3 = bondPartner[settleIatoms[i + 3]];
            if (partner3 == -1)
            {
                partner3 = settleIatoms[i + 1];
            }
            else
            {
                partner3 = -2;
            }
        }

        // Store the original masses in a temporary list
        std::vector<real> originalMasses(atoms.nr);
        for (int a = 0; a < atoms.nr; a++)
        {
            if (useFep && (atoms.atom[a].m < minMass || atoms.atom[a].mB < minMass)
                && (atoms.atom[a].m != atoms.atom[a].mB
                    || (bondPartner[a] >= 0
                        && (atoms.atom[bondPartner[a]].m != atoms.atom[bondPartner[a]].mB))))
            {
                wi->addError("Hydrogen mass repartitioning is not supported with perturbed masses");

                return;
            }

            originalMasses[a] = atoms.atom[a].m;
        }

        for (int a = 0; a < atoms.nr; a++)
        {
            if (originalMasses[a] > 0 && originalMasses[a] < minMass)
            {
                const int partner = bondPartner[a];

                if (partner >= 0)
                {
                    GMX_RELEASE_ASSERT(atoms.atom[a].m == originalMasses[a],
                                       "Here we should have unmodified masses");

                    const real massDiff = minMass - atoms.atom[a].m;

                    if (originalMasses[partner] - massDiff < minMass)
                    {
                        haveTooSmallPartnerMass = true;
                    }
                    else
                    {
                        atoms.atom[a].m  = minMass;
                        atoms.atom[a].mB = minMass;
                        atoms.atom[partner].m -= massDiff;
                        atoms.atom[partner].mB -= massDiff;

                        if (atoms.atom[partner].m < minMass)
                        {
                            haveTooSmallPartnerMass = true;
                        }
                    }
                }
                else if (partner == -1)
                {
                    numUnboundLightAtoms++;
                }
                else
                {
                    numMultipleBoundLightAtoms++;
                }
            }
        }
    }

    if (haveTooSmallPartnerMass)
    {
        wi->addError(
                "Light atoms are bound to at least one atom that has a too low mass for "
                "repartioning");
    }
    if (numUnboundLightAtoms > 0)
    {
        wi->addWarning(
                formatString("The are %d atoms that have a mass below the mass repartitioning "
                             "limit but are not bound. These masses cannot be repartitioned.",
                             numUnboundLightAtoms));
    }
    if (numMultipleBoundLightAtoms)
    {
        wi->addError(formatString(
                "The are %d atoms that have a mass below the mass repartitioning limit and have "
                "multiple bonds, whereas they should have only one bond.",
                numMultipleBoundLightAtoms));
    }
}

} // namespace gmx
