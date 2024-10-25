/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#ifndef GMX_GMXPREPROCESS_GEN_AD_H
#define GMX_GMXPREPROCESS_GEN_AD_H

#include "specbond.h"

struct t_atoms;
struct t_excls;
struct MoleculePatchDatabase;
struct InteractionsOfType;
struct PreprocessResidue;
struct DisulfideBond;

namespace gmx
{
template<typename>
class ArrayRef;
}

/*! \brief
 * Generate pairs, angles and dihedrals from .rtp settings
 *
 * \param[in,out] atoms            Global information about atoms in topology.
 * \param[in]     rtpFFDB          Residue type database from force field.
 * \param[in,out] plist            Information about listed interactions.
 * \param[in,out] excls            Pair interaction exclusions.
 * \param[in,out] globalPatches    Information about possible residue modifications.
 * \param[in]     bAllowMissing    True if missing interaction information is allowed.
 *                                 AKA allow cartoon physics
 * \param[in]     cyclicBondsIndex Information about bonds creating cyclic molecules.
 *                                 Empty if no such bonds exist.
 * \param[in]     ssbonds          Information regarding special bonds and custom
 *                                 improper dihedrals from specbond. Empty if no such
 *                                 bonds exist.
 */
void gen_pad(t_atoms*                               atoms,
             gmx::ArrayRef<const PreprocessResidue> rtpFFDB,
             gmx::ArrayRef<InteractionsOfType>      plist,
             t_excls                                excls[],
             gmx::ArrayRef<MoleculePatchDatabase>   globalPatches,
             bool                                   bAllowMissing,
             gmx::ArrayRef<const int>               cyclicBondsIndex,
             gmx::ArrayRef<const DisulfideBond>     ssbonds);

#endif
