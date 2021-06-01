/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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

namespace gmx
{
class MDLogger;
template<typename>
class ArrayRef;
} // namespace gmx

namespace gmx
{

/*! \brief Has responsibility for checking that the local topology
 * distributed across domains describes a total number of bonded
 * interactions that matches the system topology
 *
 * Because this check is not urgent, the communication that it
 * requires is done at the next opportunity, rather than requiring
 * extra communication. If the check fails, a fatal error stops
 * execution. In principle, if there was a bug, GROMACS might crash in
 * the meantime because of the wrong forces. However as a bug is
 * unlikely, we optimize by avoiding creating extra overhead from
 * communication.
 */
class LocalTopologyChecker
{
public:
    /*! \brief Constructor
     * \param[in]    mdlog            Logger
     * \param[in]    cr               Communication object
     * \param[in]    mtop             Global system topology
     * \param[in]    useUpdateGroups  Whether update groups are in use
     */
    LocalTopologyChecker(const MDLogger& mdlog, const t_commrec* cr, const gmx_mtop_t& mtop, bool useUpdateGroups);
    //! Destructor
    ~LocalTopologyChecker();
    //! Move constructor
    LocalTopologyChecker(LocalTopologyChecker&& other) noexcept;
    //! Move assignment
    LocalTopologyChecker& operator=(LocalTopologyChecker&& other) noexcept;

    /*! \brief Set that the local topology should be checked via
     * observables reduction whenever that reduction is required by
     * another module. */
    void scheduleCheckOfLocalTopology(int numBondedInteractionsToReduce);

    /*! \brief Return whether the total bonded interaction count across
     * domains should be checked in observables reduction this step. */
    bool shouldCheckNumberOfBondedInteractions() const;

    //! Return the number of bonded interactions in this domain.
    int numBondedInteractions() const;

    /*! \brief Set total bonded interaction count across domains. */
    void setNumberOfBondedInteractionsOverAllDomains(int newValue);

    /*! \brief Check whether bonded interactions are missing from the reverse topology
     * produced by domain decomposition.
     *
     * \param[in]    top_local  Local topology for the error message
     * \param[in]    x          Position vector for the error message
     * \param[in]    box        Box matrix for the error message
     */
    void checkNumberOfBondedInteractions(const gmx_localtop_t* top_local,
                                         ArrayRef<const RVec>  x,
                                         const matrix          box);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx
#endif
