/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Declares the functions for generating update groups
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_UPDATEGROUPS
#define GMX_MDLIB_UPDATEGROUPS

#include <variant>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

enum class PbcType : int;
struct gmx_mtop_t;

namespace gmx
{
class MDLogger;
class RangePartitioning;

/*! \brief Returns a vector with update groups for each moleculetype in \p mtop
 * or an error string when the criteria (see below) are not satisfied.
 *
 * An error string is returned when at least one moleculetype does not obey
 * the restrictions of update groups, e.g. more than two constraints in a row.
 *
 * Currently valid update groups are:
 * - a single atom which is not a virtual site and does not have constraints;
 * - or a group of atoms where all virtual sites are constructed from atoms
 *   within the group and at least one non-vsite atom is constrained to
 *   all other non-vsite atoms.
 *
 * To have update groups, all virtual sites should be linear 2 or 3 atom
 * constructions with coefficients >= 0 and sum of coefficients <= 1.
 *
 * This vector is generally consumed in constructing an UpdateGroups object.
 *
 * \param[in] mtop  The system topology
 */
std::variant<std::vector<RangePartitioning>, std::string> makeUpdateGroupingsPerMoleculeType(const gmx_mtop_t& mtop);

/*! \brief Returns the maximum update group radius
 *
 * \note When \p updateGroups is empty, 0 is returned.
 *
 * \param[in] mtop                           The system topology
 * \param[in] updateGroupingPerMoleculeType  List of update group, size should match the
 *                                           number of moltypes in \p mtop or be 0
 * \param[in] temperature                    The maximum reference temperature, pass -1
 *                                           when unknown or not applicable
 */
real computeMaxUpdateGroupRadius(const gmx_mtop_t&                 mtop,
                                 ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType,
                                 real                              temperature);

/*! \brief Returns whether mtop contains any constraints and/or vsites
 *
 * When we have constraints and/or vsites, it is beneficial to use
 * update groups (when possible) to allow independent update of
 * groups.*/
bool systemHasConstraintsOrVsites(const gmx_mtop_t& mtop);

/*! \libinternal
 * \brief Owns the update grouping and related data */
class UpdateGroups
{
public:
    //! Default constructor
    UpdateGroups() = default;
    //! Constructor when update groups are active
    UpdateGroups(std::vector<RangePartitioning>&& updateGroupingPerMoleculeType, real maxUpdateGroupRadius);

    bool                              useUpdateGroups() const { return useUpdateGroups_; }
    real                              maxUpdateGroupRadius() const { return maxUpdateGroupRadius_; }
    ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType() const
    {
        return updateGroupingPerMoleculeType_;
    }

private:
    //! Whether update groups are in use
    bool useUpdateGroups_ = false;
    //! The update groupings within each respective molecule type, empty when not in use
    std::vector<RangePartitioning> updateGroupingPerMoleculeType_ = {};
    //! The maximum radius of any update group, 0 when not in use
    real maxUpdateGroupRadius_ = 0.0_real;
};

/*! \brief Builder for update groups.
 *
 * Checks the conditions for using update groups, and logs a message
 * if they cannot be used, along with the reason why not.
 *
 * \p updateGroupingPerMoleculeType can not be empty (is asserted).
 *
 * If PP domain decomposition is not in use, there is no reason to use
 * update groups.
 *
 * All molecule types in the system topology must be conform to the
 * requirements, such that makeUpdateGroupingsPerMoleculeType()
 * returns a non-empty vector.
 *
 * mdrun -rerun does not support update groups (PBC corrections needed).
 *
 * When we have constraints and/or vsites, it is beneficial to use
 * update groups (when possible) to allow independent update of
 * groups. But if there are no constraints or vsites, then there is no
 * need to use update groups at all.
 *
 * To use update groups, the large domain-to-domain cutoff distance
 * should be compatible with the box size.
 */
UpdateGroups makeUpdateGroups(const gmx::MDLogger&             mdlog,
                              std::vector<RangePartitioning>&& updateGroupingPerMoleculeType,
                              real                             maxUpdateGroupRadius,
                              bool                             doRerun,
                              bool                             useDomainDecomposition,
                              bool                             systemHasConstraintsOrVsites,
                              real                             cutoffMargin);

} // namespace gmx

#endif // GMX_MDLIB_UPDATEGROUPS
