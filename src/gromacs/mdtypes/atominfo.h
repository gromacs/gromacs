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
/*! \libinternal \file
 *
 * \brief This file makes declarations used for storing bitfields
 * describing each atom so that other modules can efficiently process
 * them.
 *
 * \inlibraryapi
 */

#ifndef GMX_MDTYPES_ATOMINFO_H
#define GMX_MDTYPES_ATOMINFO_H

#include <cstdint>

#include <vector>

namespace gmx
{

/*! \brief Constants whose bit describes a property of an atom in
 * AtomInfoWithinMoleculeBlock.atomInfo.
 *
 * No bit should exceed 1 << 31, so that it fits into a 32-bit
 * integer.
 *
 * Since the tpx format support max 256 energy groups, we do the same
 * here, reserving bits 0-7 for the energy-group ID.
 */
//! \{
static constexpr int32_t sc_atomInfo_FreeEnergyPerturbation = 1 << 15;
static constexpr int32_t sc_atomInfo_HasPerturbedCharge     = 1 << 16;
static constexpr int32_t sc_atomInfo_Exclusion              = 1 << 17;
static constexpr int32_t sc_atomInfo_Constraint             = 1 << 20;
static constexpr int32_t sc_atomInfo_Settle                 = 1 << 21;
static constexpr int32_t sc_atomInfo_BondCommunication      = 1 << 22;
static constexpr int32_t sc_atomInfo_HasVdw                 = 1 << 23;
static constexpr int32_t sc_atomInfo_HasCharge              = 1 << 24;
//! \}
//! The first 8 bits are reserved for energy-group ID
static constexpr int32_t sc_atomInfo_EnergyGroupIdMask = 0b11111111;

/*! \internal
 *  \brief Contains information about each atom in a molecule block of the global topology.
 */
struct AtomInfoWithinMoleculeBlock
{
    //! Index within the system of the first atom in the molecule block
    int indexOfFirstAtomInMoleculeBlock = 0;
    //! Index within the system of the last atom in the molecule block
    int indexOfLastAtomInMoleculeBlock = 0;
    /*! \brief Atom info for each atom in the block.
     *
     * The typical case is that all atoms are identical for each
     * molecule of the block, and if so this vector has size equal to
     * the number of atoms in the molecule.
     *
     * An example of an atypical case is QM/MM, where multiple
     * molecules might be present and different molecules have
     * different atoms within any one QM group or region. Now there are
     * multiple kinds of molecules with the same connectivity, so we simply
     * write out the atom info for the entire molecule block. Then the
     * size equals the product of the number of atoms in the
     * molecule and the number of molecules.
     *
     * The vector needs to be indexed accordingly.
     */
    std::vector<int32_t> atomInfo;
};

} // namespace gmx

#endif
