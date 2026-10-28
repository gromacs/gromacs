/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Defines builder function for handling PDB CONECT records
 *
 * \ingroup topology
 */
#include "gmxpre.h"

#include "gromacs/topology/pdb_conect.h"

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

PdbConect buildPdbConect(const gmx_mtop_t& mtop)
{
    PdbConect pdbConect;

    int indexOfFirstAtomInMolecule = 0;
    // Loop over the molecule blocks
    for (const gmx_molblock_t& molblock : mtop.molblock)
    {
        const gmx_moltype_t& moltype = mtop.moltype[molblock.type];
        // Loop over each molecule in the block
        for (int mol = 0; mol < molblock.nmol; mol++, indexOfFirstAtomInMolecule += moltype.atoms.nr)
        {
            // Find each "chemical bond"
            for (const auto functionType : gmx::EnumerationWrapper<InteractionFunction>{})
            {
                if (IS_CHEMBOND(functionType))
                {
                    // Record a connection between its two atoms
                    const InteractionList& ilist = moltype.ilist[functionType];
                    for (int i = 1; i < ilist.size(); i += NRAL(functionType) + 1)
                    {
                        PdbConectRecord newValue{ indexOfFirstAtomInMolecule + ilist.iatoms[i],
                                                  indexOfFirstAtomInMolecule + ilist.iatoms[i + 1] };
                        pdbConect.insert(newValue);
                    }
                }
            }
        }
    }

    return pdbConect;
}

} // namespace gmx
