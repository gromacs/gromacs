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
 * Declarations for chain ID construction routines.
 *
 * \ingroup module_topology
 */
#ifndef GMX_TOPOLOGY_CHAINID_H
#define GMX_TOPOLOGY_CHAINID_H

struct t_atoms;
struct t_block;

namespace gmx
{

class RangePartitioning;

/*! \brief Make chain IDs in \p atoms for legacy t_block \p mols.
 *
 * Chain IDs will only be constructed if \p mols contains at least 2 chains
 * with at least 15 particles each.
 *
 * Regardless of whether new chain IDs are constructed, already present
 * chain IDs in \p atoms may be cleared (set to ' ').
 *
 * \note Chain IDs can and will only be constructed in \p atoms if \p mols was
 * created from a .tpr file, i.e. with the full system topology.
 */
void tpx_make_chain_identifiers(t_atoms* atoms, const t_block* mols);

/*! \brief Make chain IDs in \p atoms for range partitioning \p mols.
 *
 * Chain IDs will only be constructed if \p mols contains at least 2 chains
 * with at least 15 particles each.
 *
 * Regardless of whether new chain IDs are constructed, already present
 * chain IDs in \p atoms may be cleared (set to ' ').
 *
 * \note Chain IDs can and will only be constructed in \p atoms if \p mols was
 * created from a .tpr file, i.e. with the full system topology.
 */
void makeChainIdentifiersAfterTprReading(t_atoms* atoms, const RangePartitioning& mols);

} // namespace gmx

#endif
