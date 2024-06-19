/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 *
 * \brief Defines a function that makes the list of links between
 * atoms connected by bonded interactions.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "gromacs/domdec/makebondedlinks.h"

#include <cstdio>

#include <algorithm>
#include <memory>
#include <vector>

#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/options.h"
#include "gromacs/domdec/reversetopology.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"

using gmx::ArrayRef;
using gmx::DDBondedChecking;

/*! \brief Check if a link is stored in \p links to \p atom and if this is not the case, store a link */
static void check_link(std::vector<int>* links, const int atom)
{
    const auto it = find(links->begin(), links->end(), atom);
    if (it == links->end())
    {
        links->push_back(atom);
    }
}

/*! \brief Creates and return a list of bonded links for all atoms in the system */
static gmx::ListOfLists<int> genBondedLinks(const gmx_mtop_t&                               mtop,
                                            gmx::ArrayRef<gmx::AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock)
{
    gmx::ListOfLists<int> link;

    /* For each atom make a list of other atoms in the system
     * that a linked to it via bonded interactions
     * which are also stored in reverse_top.
     */

    reverse_ilist_t ril_intermol;
    if (mtop.bIntermolecularInteractions)
    {
        t_atoms atoms;

        atoms.nr   = mtop.natoms;
        atoms.atom = nullptr;

        GMX_RELEASE_ASSERT(mtop.intermolecular_ilist,
                           "We should have an ilist when intermolecular interactions are on");

        ReverseTopOptions rtOptions(DDBondedChecking::ExcludeZeroLimit);
        make_reverse_ilist(
                *mtop.intermolecular_ilist, &atoms, rtOptions, AtomLinkRule::AllAtomsInBondeds, &ril_intermol);
    }

    int indexOfFirstAtomInMolecule = 0;
    int numLinkedAtoms             = 0;
    for (size_t mb = 0; mb < mtop.molblock.size(); mb++)
    {
        const gmx_molblock_t& molb = mtop.molblock[mb];
        if (molb.nmol == 0)
        {
            continue;
        }
        const gmx_moltype_t& molt = mtop.moltype[molb.type];
        /* Make a reverse ilist in which the interactions are linked
         * to all atoms, not only the first atom as in gmx_reverse_top.
         * The constraints are discarded here.
         */
        ReverseTopOptions rtOptions(DDBondedChecking::ExcludeZeroLimit);
        reverse_ilist_t   ril;
        make_reverse_ilist(molt.ilist, &molt.atoms, rtOptions, AtomLinkRule::AllAtomsInBondeds, &ril);

        gmx::AtomInfoWithinMoleculeBlock* atomInfoOfMoleculeBlock = &atomInfoForEachMoleculeBlock[mb];

        std::vector<int> linksForOneAtom;
        int              mol = 0;
        for (mol = 0; mol < (mtop.bIntermolecularInteractions ? molb.nmol : 1); mol++)
        {
            for (int a = 0; a < molt.atoms.nr; a++)
            {
                const int atomIndex = indexOfFirstAtomInMolecule + a;
                linksForOneAtom.clear();
                int i = ril.index[a];
                while (i < ril.index[a + 1])
                {
                    int ftype = ril.il[i++];
                    int nral  = NRAL(ftype);
                    /* Skip the ifunc index */
                    i++;
                    for (int j = 0; j < nral; j++)
                    {
                        int aj = ril.il[i + j];
                        if (aj != a)
                        {
                            check_link(&linksForOneAtom, indexOfFirstAtomInMolecule + aj);
                        }
                    }
                    i += nral_rt(ftype);
                }

                if (mtop.bIntermolecularInteractions)
                {
                    int i = ril_intermol.index[atomIndex];
                    while (i < ril_intermol.index[atomIndex + 1])
                    {
                        int ftype = ril_intermol.il[i++];
                        int nral  = NRAL(ftype);
                        /* Skip the ifunc index */
                        i++;
                        for (int j = 0; j < nral; j++)
                        {
                            /* Here we assume we have no charge groups;
                             * this has been checked above.
                             */
                            int aj = ril_intermol.il[i + j];
                            check_link(&linksForOneAtom, aj);
                        }
                        i += nral_rt(ftype);
                    }
                }
                if (!linksForOneAtom.empty())
                {
                    atomInfoOfMoleculeBlock->atomInfo[a] |= gmx::sc_atomInfo_BondCommunication;
                    numLinkedAtoms++;
                }
                // Add the links for the current atom to the total list
                link.pushBack(linksForOneAtom);
            }

            indexOfFirstAtomInMolecule += molt.atoms.nr;

            GMX_ASSERT(link.size() == size_t(indexOfFirstAtomInMolecule),
                       "We should have one link entry for each atom");
        }
        int nlink_mol = link.listRangesView()[indexOfFirstAtomInMolecule]
                        - link.listRangesView()[indexOfFirstAtomInMolecule - molt.atoms.nr];

        if (debug)
        {
            fprintf(debug,
                    "molecule type '%s' %d atoms has %d atom links through bonded interac.\n",
                    *molt.name,
                    molt.atoms.nr,
                    nlink_mol);
        }

        if (molb.nmol > mol)
        {
            /* Copy the data for the rest of the molecules in this block */
            for (; mol < molb.nmol; mol++)
            {
                for (int a = 0; a < molt.atoms.nr; a++)
                {
                    const int atomIndex = indexOfFirstAtomInMolecule + a;
                    link.pushBackListOfSize(link[atomIndex - molt.atoms.nr].size());
                    const gmx::ArrayRef<const int> linksForAtomInPreviousMolecule =
                            link[atomIndex - molt.atoms.nr];
                    gmx::ArrayRef<int> linksForAtom = link.back();
                    std::transform(linksForAtomInPreviousMolecule.begin(),
                                   linksForAtomInPreviousMolecule.end(),
                                   linksForAtom.begin(),
                                   [&molt](const auto a) { return a + molt.atoms.nr; });
                    if (!linksForAtom.empty()
                        && atomIndex - atomInfoOfMoleculeBlock->indexOfFirstAtomInMoleculeBlock
                                   < gmx::ssize(atomInfoOfMoleculeBlock->atomInfo))
                    {
                        atomInfoOfMoleculeBlock->atomInfo[atomIndex - atomInfoOfMoleculeBlock->indexOfFirstAtomInMoleculeBlock] |=
                                gmx::sc_atomInfo_BondCommunication;
                        numLinkedAtoms++;
                    }
                }
                indexOfFirstAtomInMolecule += molt.atoms.nr;
            }
        }
    }

    if (debug)
    {
        fprintf(debug, "Of the %d atoms %d are linked via bonded interactions\n", mtop.natoms, numLinkedAtoms);
    }

    return link;
}

void makeBondedLinks(gmx_domdec_t*                                   dd,
                     const gmx_mtop_t&                               mtop,
                     gmx::ArrayRef<gmx::AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock)
{
    if (dd->comm->systemInfo.filterBondedCommunication)
    {
        dd->comm->bondedLinks  = std::make_unique<gmx::ListOfLists<int>>();
        *dd->comm->bondedLinks = genBondedLinks(mtop, atomInfoForEachMoleculeBlock);
    }
}
