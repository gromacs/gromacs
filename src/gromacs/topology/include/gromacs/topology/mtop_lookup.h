/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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

/*! Class for quickly looking up properties of atoms based on global atom index
 *
 * This class stores the index of the molecule block for the last lookup, which
 * enables faster lookup for series of indices belonging to the same molecule block.
 */
class MTopLookUp
{
public:
    /*! \brief Constructor
     *
     * \param[in] mtop  A reference to the molecule topology
     */
    MTopLookUp(const gmx_mtop_t& mtop) : mtop_(mtop) {}

    //! Struct for returning the molecular topology indices of a global atom index
    struct MolblockAtomIndex
    {
        //! The molecule block index in \p mtop
        int molBlock;
        //! The index of the molecule in the block
        int molIndex;
        //! The atom index in the molecule
        int atomIndex;
    };

    /*! \brief Look up the molecule block and other indices of a global atom index
     *
     * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
     *
     * \param[in] globalAtomIndex  The global atom index to look up
     *
     * \returns the molecular topology indices
     */
    MolblockAtomIndex getMolblockAtomIndex(const int globalAtomIndex)
    {
        GMX_ASSERT(globalAtomIndex >= 0, "The atom index to look up should not be negative");
        GMX_ASSERT(globalAtomIndex < mtop_.natoms,
                   "The atom index to look up should be within range");
        GMX_ASSERT(!mtop_.moleculeBlockIndices.empty(),
                   "The moleculeBlockIndices should not be empty");

        /* Search the molecule block index using bisection */
        int molBlock0 = -1;
        int molBlock1 = mtop_.molblock.size();

        int globalAtomStart = 0;
        while (true)
        {
            globalAtomStart = mtop_.moleculeBlockIndices[molBlock_].globalAtomStart;
            if (globalAtomIndex < globalAtomStart)
            {
                molBlock1 = molBlock_;
            }
            else if (globalAtomIndex >= mtop_.moleculeBlockIndices[molBlock_].globalAtomEnd)
            {
                molBlock0 = molBlock_;
            }
            else
            {
                break;
            }
            molBlock_ = ((molBlock0 + molBlock1 + 1) >> 1);
        }

        int molIndex = (globalAtomIndex - globalAtomStart)
                       / mtop_.moleculeBlockIndices[molBlock_].numAtomsPerMolecule;

        return { molBlock_,
                 molIndex,
                 globalAtomIndex - globalAtomStart
                         - molIndex * mtop_.moleculeBlockIndices[molBlock_].numAtomsPerMolecule };
    }

    /*! \brief Returns the global molecule index of a global atom index
     *
     * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
     *
     * \param[in] globalAtomIndex  The global atom index to look up
     */
    int getMoleculeIndex(const int globalAtomIndex)
    {
        const MolblockAtomIndex mbai = getMolblockAtomIndex(globalAtomIndex);

        return mtop_.moleculeBlockIndices[mbai.molBlock].moleculeIndexStart + mbai.molIndex;
    }

    /*! \brief Returns the atom data for an atom based on global atom index
     *
     * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
     *
     * \param[in] globalAtomIndex  The global atom index to look up
     */
    const t_atom& getAtomParameters(const int globalAtomIndex)
    {
        const MolblockAtomIndex mbai = getMolblockAtomIndex(globalAtomIndex);

        const gmx_moltype_t& moltype = mtop_.moltype[mtop_.molblock[mbai.molBlock].type];
        return moltype.atoms.atom[mbai.atomIndex];
    }

    //! Struct for returning atom and residue name and index
    struct AtomAndResidueNameAndNumber
    {
        //! The name of the atom
        const char* atomName;
        //! The residue number of the atom in the molecule
        int residueNumber;
        //! The residue name of the atom
        const char* residueName;
    };

    /*! \brief Look up the atom and residue name and residue number and index of a global atom index
     *
     * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
     * Note that this function does a (somewhat expensive) lookup. If you want
     * to look up data sequentially for all atoms in a molecule or the system,
     * use one of the mtop loop functionalities.
     *
     * \param[in] globalAtomIndex  The global atom index to look up
     *
     * \returns a struct with atom and residue information
     */
    AtomAndResidueNameAndNumber getAtomAndResidueNameAndNumber(const int globalAtomIndex)
    {
        const MolblockAtomIndex mbai = getMolblockAtomIndex(globalAtomIndex);

        const t_atoms&              atoms = mtop_.moltype[mtop_.molblock[mbai.molBlock].type].atoms;
        const MoleculeBlockIndices& indices  = mtop_.moleculeBlockIndices[mbai.molBlock];
        const int                   resIndex = atoms.atom[mbai.atomIndex].resind;

        AtomAndResidueNameAndNumber aarnan;

        aarnan.atomName = *(atoms.atomname[mbai.atomIndex]);
        if (atoms.nres > mtop_.maxResiduesPerMoleculeToTriggerRenumber())
        {
            aarnan.residueNumber = atoms.resinfo[resIndex].nr;
        }
        else
        {
            /* Single residue molecule, keep counting */
            aarnan.residueNumber = indices.residueNumberStart + mbai.molIndex * atoms.nres + resIndex;
        }
        aarnan.residueName = *(atoms.resinfo[resIndex].name);

        return aarnan;
    }

    /*! \brief Look up global residue index for a global atom index
     *
     * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
     * Note that this function does a (somewhat expensive) lookup. If you want
     * to look up data sequentially for all atoms in a molecule or the system,
     * use one of the mtop loop functionalities.
     *
     * \param[in] globalAtomIndex  The global atom index to look up
     */
    int getGlobalResidueIndex(const int globalAtomIndex)
    {
        const MolblockAtomIndex mbai = getMolblockAtomIndex(globalAtomIndex);

        const t_atoms&              atoms = mtop_.moltype[mtop_.molblock[mbai.molBlock].type].atoms;
        const MoleculeBlockIndices& indices = mtop_.moleculeBlockIndices[mbai.molBlock];

        return indices.globalResidueStart + mbai.molIndex * atoms.nres + atoms.atom[mbai.atomIndex].resind;
    }

    /*! \brief Returns residue information for an atom based on global atom index
     *
     * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
     *
     * \param[in] globalAtomIndex  The global atom index to look up
     */
    const t_resinfo& getResidueInfo(const int globalAtomIndex)
    {
        const MolblockAtomIndex mbai = getMolblockAtomIndex(globalAtomIndex);

        const gmx_moltype_t& moltype = mtop_.moltype[mtop_.molblock[mbai.molBlock].type];
        const int            resind  = moltype.atoms.atom[mbai.atomIndex].resind;

        return moltype.atoms.resinfo[resind];
    }

    /*! \brief Returns PDB information for an atom based on global atom index
     *
     * The atom index has to be in range: 0 <= \p globalAtomIndex < \p mtop->natoms.
     *
     * \param[in] globalAtomIndex  The global atom index to look up
     */
    const t_pdbinfo& getAtomPdbInfo(const int globalAtomIndex)
    {
        const MolblockAtomIndex mbai = getMolblockAtomIndex(globalAtomIndex);

        const gmx_moltype_t& moltype = mtop_.moltype[mtop_.molblock[mbai.molBlock].type];
        GMX_ASSERT(moltype.atoms.havePdbInfo, "PDB information not present when requested");

        return moltype.atoms.pdbinfo[mbai.atomIndex];
    }

    //! Returns a reference to the molecular topology
    const gmx_mtop_t& mtop() const { return mtop_; }

private:
    //! A reference to the molecular topology
    const gmx_mtop_t& mtop_;
    //! The last molBlock index that was found, stored for faster look up of the next index
    int molBlock_ = 0;
};

#endif
