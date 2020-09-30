/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
#ifndef GMX_MULTIPLETIMESTEPPING_H
#define GMX_MULTIPLETIMESTEPPING_H

#include <bitset>

#include "gromacs/utility/enumerationhelpers.h"

struct t_inputrec;

namespace gmx
{

//! Force group available for selection for multiple time step integration
enum class MtsForceGroups : int
{
    LongrangeNonbonded, //!< PME-mesh or Ewald for electrostatics and/or LJ
    Nonbonded,          //!< Non-bonded pair interactions
    Pair,               //!< Bonded pair interactions
    Dihedral,           //!< Dihedrals, including cmap (not restraints)
    Angle,              //! Bonded angle potentials (not restraints)
    Count               //! The number of groups above
};

static const gmx::EnumerationArray<MtsForceGroups, std::string> mtsForceGroupNames = {
    "longrange-nonbonded", "nonbonded", "pair", "dihedral", "angle"
};

//! Setting for a single level for multiple time step integration
struct MtsLevel
{
    //! The force group selection for this level;
    std::bitset<static_cast<int>(MtsForceGroups::Count)> forceGroups;
    //! The factor between the base, fastest, time step and the time step for this level
    int stepFactor;
};

/*! \brief Returns the interval in steps at which the non-bonded pair forces are calculated
 *
 * Note: returns 1 when multiple time-stepping is not activated.
 */
int nonbondedMtsFactor(const t_inputrec& ir);

//! (Release) Asserts that all multiple time-stepping requirements on \p ir are fulfilled
void assertMtsRequirements(const t_inputrec& ir);

} // namespace gmx

#endif /* GMX_MULTIPLETIMESTEPPING_H */
