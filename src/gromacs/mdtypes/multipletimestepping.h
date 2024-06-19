/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
#ifndef GMX_MULTIPLETIMESTEPPING_H
#define GMX_MULTIPLETIMESTEPPING_H

#include <bitset>
#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/gmxassert.h"

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
    Angle,              //!< Bonded angle potentials (not restraints)
    Pull,               //!< COM pulling
    Awh,                //!< Accelerated weight histogram method
    Count               //!< The number of groups above
};

//! Names for the MTS force groups
static const gmx::EnumerationArray<MtsForceGroups, std::string> mtsForceGroupNames = {
    "longrange-nonbonded", "nonbonded", "pair", "dihedral", "angle", "pull", "awh"
};

//! Setting for a single level for multiple time step integration
struct MtsLevel
{
    //! The force group selection for this level;
    std::bitset<static_cast<int>(MtsForceGroups::Count)> forceGroups;
    //! The factor between the base, fastest, time step and the time step for this level
    int stepFactor;
};

/*! \brief Returns the MTS level at which a force group is to be computed
 *
 * \param[in] mtsLevels  List of force groups for each MTS level, can be empty without MTS
 * \param[in] mtsForceGroup  The force group to query the MTS level for
 */
static inline int forceGroupMtsLevel(ArrayRef<const MtsLevel> mtsLevels, const MtsForceGroups mtsForceGroup)
{
    GMX_ASSERT(mtsLevels.empty() || mtsLevels.size() == 2, "Only 0 or 2 MTS levels are supported");

    return (mtsLevels.empty() || mtsLevels[0].forceGroups[static_cast<int>(mtsForceGroup)]) ? 0 : 1;
};

/*! \brief Returns the interval in steps at which the non-bonded pair forces are calculated
 *
 * Note: returns 1 when multiple time-stepping is not activated.
 */
int nonbondedMtsFactor(const t_inputrec& ir);

//! Struct for passing the MTS mdp options to setupMtsLevels()
struct GromppMtsOpts
{
    //! The number of MTS levels
    int numLevels = 0;
    //! The names of the force groups assigned by the user to level 2, internal index 1
    std::string level2Forces;
    //! The step factor assigned by the user to level 2, internal index 1
    int level2Factor = 0;
};

/*! \brief Sets up and returns the MTS levels and checks requirements of MTS
 *
 * Appends errors about allowed input values ir to errorMessages, when not nullptr.
 *
 * \param[in]     mtsOpts        Options for setting the MTS levels
 * \param[in,out] errorMessages  List of error messages, can be nullptr
 */
std::vector<MtsLevel> setupMtsLevels(const GromppMtsOpts& mtsOpts, std::vector<std::string>* errorMessages);

/*! \brief Returns whether we use MTS and the MTS setup is internally valid
 *
 * Note that setupMtsLevels would have returned at least one error message
 * when this function returns false
 */
bool haveValidMtsSetup(const t_inputrec& ir);

/*! \brief Checks whether the MTS requirements on other algorithms and output frequencies are met
 *
 * Note: exits with an assertion failure when
 * ir.useMts == true && haveValidMtsSetup(ir) == false
 *
 * \param[in] ir  Complete input record
 * \returns list of error messages, empty when all MTS requirements are met
 */
std::vector<std::string> checkMtsRequirements(const t_inputrec& ir);

} // namespace gmx

#endif /* GMX_MULTIPLETIMESTEPPING_H */
