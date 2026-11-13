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
/*! \file
 * \brief
 * Definitions for chain ID construction routines.
 *
 * \ingroup module_topology
 */

#include "gmxpre.h"

#include "gromacs/topology/chainid.h"

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"

namespace gmx
{
namespace
{

//! Constructs plausible chain IDs for multi-molecule systems, e.g. when read from .tpr files
class ChainIdFiller
{
public:
    //! Fill in the chain ID for the indicated atom range, which might be a molecule.
    void fill(t_atoms* atoms, int startAtomIndex, int endAtomIndex);
    //! If only one chain was found, we don't add a chain ID.
    void clearIfNeeded(t_atoms* atoms) const;

private:
    //! Minimum size for a chain worth giving an ID
    static constexpr int s_chainMinAtoms = 15;

    //! The number of the next chain that will be assigned.
    int nextChainNumber_ = 0;
    //! The chain ID of the next chain that will be assigned.
    char nextChainId_ = 'A';
    //! Whether the set of chain IDs (ie. upper- and lower-case letters and single digits) is exhausted.
    bool outOfIds_ = false;
};

void ChainIdFiller::fill(t_atoms* atoms, const int startAtomIndex, const int endAtomIndex)
{
    // We always assign a new chain number, but only assign a chain id
    // characters for larger molecules.
    char chainIdToAssign;
    if (endAtomIndex - startAtomIndex >= s_chainMinAtoms && !outOfIds_)
    {
        /* Set the chain id for the output */
        chainIdToAssign = nextChainId_;
        /* Here we allow for the max possible 2*26+10=62 chain ids */
        if (nextChainId_ == 'Z')
        {
            nextChainId_ = 'a';
        }
        else if (nextChainId_ == 'z')
        {
            nextChainId_ = '0';
        }
        else if (nextChainId_ == '9')
        {
            outOfIds_ = true;
        }
        else
        {
            nextChainId_++;
        }
    }
    else
    {
        chainIdToAssign = ' ';
    }
    for (int a = startAtomIndex; a < endAtomIndex; a++)
    {
        atoms->resinfo[atoms->atom[a].resind].chainnum = nextChainNumber_;
        atoms->resinfo[atoms->atom[a].resind].chainid  = chainIdToAssign;
    }
    nextChainNumber_++;
}

void ChainIdFiller::clearIfNeeded(t_atoms* atoms) const
{
    /* Blank out the chain id if there was only one chain */
    if (nextChainId_ == 'B')
    {
        for (int r = 0; r < atoms->nres; r++)
        {
            atoms->resinfo[r].chainid = ' ';
        }
    }
}

}; // namespace

void makeChainIdentifiersAfterTprReading(t_atoms* atoms, const RangePartitioning& mols)
{
    ChainIdFiller filler;
    for (auto m = 0; m != mols.numBlocks(); ++m)
    {
        filler.fill(atoms, mols.block(m).begin(), mols.block(m).end());
    }
    filler.clearIfNeeded(atoms);
}

void tpx_make_chain_identifiers(t_atoms* atoms, const t_block* mols)
{
    ChainIdFiller filler;
    for (int m = 0; m < mols->nr; m++)
    {
        filler.fill(atoms, mols->index[m], mols->index[m + 1]);
    }
    filler.clearIfNeeded(atoms);
}

} // namespace gmx
