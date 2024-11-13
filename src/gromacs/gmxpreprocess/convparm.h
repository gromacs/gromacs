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

#ifndef GMX_GMXPREPROCESS_CONVPARM_H
#define GMX_GMXPREPROCESS_CONVPARM_H

#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct MoleculeInformation;
struct InteractionsOfType;
enum class CombinationRule : int;

namespace gmx
{
template<typename>
class ArrayRef;
}

//! Whether interactions of type \c ftype should have parameters converted
bool shouldConvertInteractionType(int ftype);

/*! \brief Convert the raw arrays of parameters read from the topology
 * file to those which match the interaction type and fill into \c mtop.
 *
 * This can involve pre-computation e.g. for G96 interactions.
 * Interaction parameter sets that are all zero are not added to \c mtop
 *
 * \param[in]  atnr    Number of atoms
 * \param[in]  nbtypes Non-bonded interaction types
 * \param[in]  mi      Molecule types read from topology file containing bonded-style interactions
 * \param[in]  intermolecular_interactions  Molecule type for inter-molecular interactions
 * \param[in]  comb    The combination rule
 * \param[in]  reppow  The repulsion power
 * \param[in]  fudgeQQ The 1-4 scaling factor
 * \param[out] mtop    The molecular topology containing the converted parameters
 */
void convertInteractionsOfType(int                                      atnr,
                               gmx::ArrayRef<const InteractionsOfType>  nbtypes,
                               gmx::ArrayRef<const MoleculeInformation> mi,
                               const MoleculeInformation*               intermolecular_interactions,
                               CombinationRule                          comb,
                               double                                   reppow,
                               real                                     fudgeQQ,
                               gmx_mtop_t*                              mtop);

#endif
