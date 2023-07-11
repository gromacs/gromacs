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
 * \brief This file declares functionality for checking whether
 * local topologies describe all bonded interactions.
 *
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_LOCALTOPOLOGYCHECKER_H
#define GMX_DOMDEC_LOCALTOPOLOGYCHECKER_H

#include <memory>

#include "gromacs/math/vectypes.h"

struct gmx_localtop_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_inputrec;
class t_state;

namespace gmx
{
enum class DDBondedChecking : bool;
class MDLogger;
class ObservablesReducerBuilder;
} // namespace gmx

namespace gmx
{

/*! \libinternal
 * \brief Has responsibility for checking that the local topology distributed
 * across domains describes a total number of bonded interactions that matches
 * the system topology
 *
 * This uses the ObservablesReducer framework to check that the count
 * of bonded interactions in the local topology made for each domain
 * sums to the expected value. Because this check is not urgent, the
 * communication that it requires is done at the next opportunity,
 * rather than requiring extra communication. If the check fails, a
 * fatal error stops execution. In principle, if there was a bug,
 * GROMACS might crash in the meantime because of the wrong
 * forces. However as a bug is unlikely we optimize by avoiding
 * creating extra overhead from communication.
 */
class LocalTopologyChecker
{
public:
    /*! \brief Constructor
     * \param[in]    mdlog            Logger
     * \param[in]    cr               Communication object
     * \param[in]    mtop             Global system topology
     * \param[in]    ddBondedChecking Tells for which bonded interactions presence should be checked
     * \param[in]    localTopology    The local topology
     * \param[in]    localState       The local state
     * \param[in]    useUpdateGroups  Whether update groups are in use
     * \param[in]    observablesReducerBuilder  Handle to builder for ObservablesReducer
     */
    LocalTopologyChecker(const MDLogger&            mdlog,
                         const t_commrec*           cr,
                         const gmx_mtop_t&          mtop,
                         DDBondedChecking           ddBondedChecking,
                         const gmx_localtop_t&      localTopology,
                         const t_state&             localState,
                         bool                       useUpdateGroups,
                         ObservablesReducerBuilder* observablesReducerBuilder);
    //! Destructor
    ~LocalTopologyChecker();
    //! Move constructor
    LocalTopologyChecker(LocalTopologyChecker&& other) noexcept;
    //! Move assignment
    LocalTopologyChecker& operator=(LocalTopologyChecker&& other) noexcept;

    /*! \brief Set that the local topology should be checked via
     * observables reduction whenever that reduction is required by
     * another module. In case of a single domain a direct assertion
     * is performed instead.
     */
    void scheduleCheckOfLocalTopology(int numBondedInteractionsToReduce);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx
#endif
