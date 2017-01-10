/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 *
 * \brief This file contains inline functions to look up atom information
 * using the global atom index.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_mtop
 */

#ifndef GMX_TOPOLOGY_MTOP_LOOKUP_H
#define GMX_TOPOLOGY_MTOP_LOOKUP_H

#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

struct t_atom;

/*! \brief Look up the molecule block and other indices of a global atom index
 *
 * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
 * The input value of moleculeBlock should be in range. Use 0 as starting value.
 * For subsequent calls to this function, e.g. in a loop, pass in the previously
 * returned value for best performance. Atoms in a group tend to be in the same
 * molecule(block), so this minimizes the search time.
 *
 * \param[in]     mtop                 The molecule topology
 * \param[in]     globalAtomIndex      The global atom index to look up
 * \param[in,out] moleculeBlock        The molecule block index in \p mtop
 * \param[out]    moleculeIndex        The index of the molecule in the block, can be NULL
 * \param[out]    atomIndexInMolecule  The atom index in the molecule, can be NULL
 */
static inline void
mtopGetMolblockIndex(const gmx_mtop_t *mtop,
                     int               globalAtomIndex,
                     int              *moleculeBlock,
                     int              *moleculeIndex,
                     int              *atomIndexInMolecule)
{
    GMX_ASSERT(globalAtomIndex >= 0 && globalAtomIndex < mtop->natoms, "The atom index to look up should be within range");
    GMX_ASSERT(moleculeBlock != nullptr, "molBlock can not be NULL");
    GMX_ASSERT(*moleculeBlock >= 0 && *moleculeBlock < mtop->nmolblock, "The starting molecule block index for the search should be within range");

    /* Search the molecue block index using bisection */
    int molBlock0 = -1;
    int molBlock1 = mtop->nmolblock;

    int globalAtomStart;
    while (TRUE)
    {
        globalAtomStart = mtop->molblock[*moleculeBlock].globalAtomStart;
        if (globalAtomIndex < globalAtomStart)
        {
            molBlock1 = *moleculeBlock;
        }
        else if (globalAtomIndex >= mtop->molblock[*moleculeBlock].globalAtomEnd)
        {
            molBlock0 = *moleculeBlock;
        }
        else
        {
            break;
        }
        *moleculeBlock = ((molBlock0 + molBlock1 + 1) >> 1);
    }

    int molIndex = (globalAtomIndex - globalAtomStart) / mtop->molblock[*moleculeBlock].natoms_mol;
    if (moleculeIndex != nullptr)
    {
        *moleculeIndex = molIndex;
    }
    if (atomIndexInMolecule != nullptr)
    {
        *atomIndexInMolecule = globalAtomIndex - globalAtomStart - molIndex*mtop->molblock[*moleculeBlock].natoms_mol;
    }
}

/*! \brief Returns the atom data for an atom based on global atom index
 *
 * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
 * The input value of moleculeBlock should be in range. Use 0 as starting value.
 * For subsequent calls to this function, e.g. in a loop, pass in the previously
 * returned value for best performance. Atoms in a group tend to be in the same
 * molecule(block), so this minimizes the search time.
 *
 * \param[in]     mtop                 The molecule topology
 * \param[in]     globalAtomIndex      The global atom index to look up
 * \param[in,out] moleculeBlock        The molecule block index in \p mtop
 */
static inline const t_atom &
mtopGetAtomParameters(const gmx_mtop_t  *mtop,
                      int                globalAtomIndex,
                      int               *moleculeBlock)
{
    int atomIndexInMolecule;
    mtopGetMolblockIndex(mtop, globalAtomIndex, moleculeBlock,
                         nullptr, &atomIndexInMolecule);
    const gmx_moltype_t &moltype = mtop->moltype[mtop->molblock[*moleculeBlock].type];
    return moltype.atoms.atom[atomIndexInMolecule];
}

/*! \brief Returns the mass of an atom based on global atom index
 *
 * Returns that A-state mass of the atom with global index \p globalAtomIndex.
 * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
 * The input value of moleculeBlock should be in range. Use 0 as starting value.
 * For subsequent calls to this function, e.g. in a loop, pass in the previously
 * returned value for best performance. Atoms in a group tend to be in the same
 * molecule(block), so this minimizes the search time.
 *
 * \param[in]     mtop                 The molecule topology
 * \param[in]     globalAtomIndex      The global atom index to look up
 * \param[in,out] moleculeBlock        The molecule block index in \p mtop
 */
static inline real
mtopGetAtomMass(const gmx_mtop_t  *mtop,
                int                globalAtomIndex,
                int               *moleculeBlock)
{
    const t_atom &atom = mtopGetAtomParameters(mtop, globalAtomIndex, moleculeBlock);
    return atom.m;
}

/*! \brief Look up the atom and residue name and residue number and index of a global atom index
 *
 * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
 * The input value of moleculeBlock should be in range. Use 0 as starting value.
 * For subsequent calls to this function, e.g. in a loop, pass in the previously
 * returned value for best performance. Atoms in a group tend to be in the same
 * molecule(block), so this minimizes the search time.
 * Note that this function does a (somewhat expensive) lookup. If you want
 * to look up data sequentially for all atoms in a molecule or the system,
 * use one of the mtop loop functionalities.
 *
 * \param[in]     mtop                The molecule topology
 * \param[in]     globalAtomIndex     The global atom index to look up
 * \param[in,out] moleculeBlock       The molecule block index in \p mtop
 * \param[out]    atomName            The atom name, input can be NULL
 * \param[out]    residueNumber       The residue number, input can be NULL
 * \param[out]    residueName         The residue name, input can be NULL
 * \param[out]    globalResidueIndex  The gobal residue index, input can be NULL
 */
static inline void
mtopGetAtomAndResidueName(const gmx_mtop_t  *mtop,
                          int                globalAtomIndex,
                          int               *moleculeBlock,
                          const char       **atomName,
                          int               *residueNumber,
                          const char       **residueName,
                          int               *globalResidueIndex)
{
    int moleculeIndex;
    int atomIndexInMolecule;
    mtopGetMolblockIndex(mtop, globalAtomIndex, moleculeBlock,
                         &moleculeIndex, &atomIndexInMolecule);

    const gmx_molblock_t &molb  = mtop->molblock[*moleculeBlock];
    const t_atoms        &atoms = mtop->moltype[molb.type].atoms;
    if (atomName != nullptr)
    {
        *atomName = *(atoms.atomname[atomIndexInMolecule]);
    }
    if (residueNumber != nullptr)
    {
        if (atoms.nres > mtop->maxres_renum)
        {
            *residueNumber = atoms.resinfo[atoms.atom[atomIndexInMolecule].resind].nr;
        }
        else
        {
            /* Single residue molecule, keep counting */
            *residueNumber = molb.residueNumberStart + moleculeIndex*atoms.nres + atoms.atom[atomIndexInMolecule].resind;
        }
    }
    if (residueName != nullptr)
    {
        *residueName = *(atoms.resinfo[atoms.atom[atomIndexInMolecule].resind].name);
    }
    if (globalResidueIndex != nullptr)
    {
        *globalResidueIndex = molb.globalResidueStart + moleculeIndex*atoms.nres + atoms.atom[atomIndexInMolecule].resind;
    }
}

/*! \brief Returns residue information for an atom based on global atom index
 *
 * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
 * The input value of moleculeBlock should be in range. Use 0 as starting value.
 * For subsequent calls to this function, e.g. in a loop, pass in the previously
 * returned value for best performance. Atoms in a group tend to be in the same
 * molecule(block), so this minimizes the search time.
 *
 * \param[in]     mtop                 The molecule topology
 * \param[in]     globalAtomIndex      The global atom index to look up
 * \param[in,out] moleculeBlock        The molecule block index in \p mtop
 */
static inline const t_resinfo &
mtopGetResidueInfo(const gmx_mtop_t  *mtop,
                   int                globalAtomIndex,
                   int               *moleculeBlock)
{
    int atomIndexInMolecule;
    mtopGetMolblockIndex(mtop, globalAtomIndex, moleculeBlock,
                         nullptr, &atomIndexInMolecule);
    const gmx_moltype_t &moltype = mtop->moltype[mtop->molblock[*moleculeBlock].type];
    const int            resind  = moltype.atoms.atom[atomIndexInMolecule].resind;
    return moltype.atoms.resinfo[resind];
}

/*! \brief Returns PDB information for an atom based on global atom index
 *
 * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
 * The input value of moleculeBlock should be in range. Use 0 as starting value.
 * For subsequent calls to this function, e.g. in a loop, pass in the previously
 * returned value for best performance. Atoms in a group tend to be in the same
 * molecule(block), so this minimizes the search time.
 *
 * \param[in]     mtop                 The molecule topology
 * \param[in]     globalAtomIndex      The global atom index to look up
 * \param[in,out] moleculeBlock        The molecule block index in \p mtop
 */
static inline const t_pdbinfo &
mtopGetAtomPdbInfo(const gmx_mtop_t  *mtop,
                   int                globalAtomIndex,
                   int               *moleculeBlock)
{
    int atomIndexInMolecule;
    mtopGetMolblockIndex(mtop, globalAtomIndex, moleculeBlock,
                         nullptr, &atomIndexInMolecule);
    const gmx_moltype_t &moltype = mtop->moltype[mtop->molblock[*moleculeBlock].type];
    GMX_ASSERT(moltype.atoms.havePdbInfo, "PDB information not present when requested");
    return moltype.atoms.pdbinfo[atomIndexInMolecule];
}

#endif
