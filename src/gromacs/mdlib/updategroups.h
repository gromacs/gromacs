/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief Declares the functions for generating update groups
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_UPDATEGROUPS
#define GMX_MDLIB_UPDATEGROUPS

#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct t_inputrec;

namespace gmx
{
class RangePartitioning;

/*! \brief Returns a vector with update groups for each moleculetype in \p mtop
 * or an empty vector when the criteria (see below) are not satisfied.
 *
 * An empty vector is returned when at least one moleculetype does not obey
 * the restrictions of update groups, e.g. more than two constraints in a row.
 *
 * Currently valid update groups are:
 * a single atom which is not a virtual site and does not have constraints;
 * or a group of atoms where all virtual sites are constructed from atoms
 * within the group and at least one non-vsite atom is constrained to
 * all other non-vsite atoms.
 * To have update groups, all virtual sites should be linear 2 or 3 atom
 * constructions with coefficients >= 0 and sum of coefficients <= 1.
 *
 * \param[in] mtop  The system topology
 */
std::vector<RangePartitioning> makeUpdateGroups(const gmx_mtop_t& mtop);

/*! \brief Returns the maximum update group radius
 *
 * \note When \p updateGroups is empty, 0 is returned.
 *
 * \param[in] mtop          The system topology
 * \param[in] updateGroups  List of update group, size should match the number of moltypes in \p mtop or be 0
 * \param[in] temperature   The maximum reference temperature, pass -1 when unknown or not applicable
 */
real computeMaxUpdateGroupRadius(const gmx_mtop_t&                      mtop,
                                 gmx::ArrayRef<const RangePartitioning> updateGroups,
                                 real                                   temperature);

} // namespace gmx

#endif // GMX_MDLIB_UPDATEGROUPS
