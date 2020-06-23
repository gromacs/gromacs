/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 The GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 * \brief Declares the VirtualSitesHandler class and vsite standalone functions
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 * \inlibraryapi
 */

#ifndef GMX_MDLIB_VSITE_H
#define GMX_MDLIB_VSITE_H

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

struct gmx_domdec_t;
struct gmx_mtop_t;
struct t_commrec;
struct InteractionList;
struct t_mdatoms;
struct t_nrnb;
struct gmx_wallcycle;
enum class PbcType : int;

namespace gmx
{
class RangePartitioning;

/*! \brief The start value of the vsite indices in the ftype enum
 *
 * The validity of the start and end values is checked in makeVirtualSitesHandler().
 * This is used to avoid loops over all ftypes just to get the vsite entries.
 * (We should replace the fixed ilist array by only the used entries.)
 */
static constexpr int c_ftypeVsiteStart = F_VSITE1;
//! The start and end value of the vsite indices in the ftype enum
static constexpr int c_ftypeVsiteEnd = F_VSITEN + 1;

//! Type for storing PBC atom information for all vsite types in the system
typedef std::array<std::vector<int>, c_ftypeVsiteEnd - c_ftypeVsiteStart> VsitePbc;

/*! \libinternal
 * \brief Class that handles construction of vsites and spreading of vsite forces
 */
class VirtualSitesHandler
{
public:
    //! Constructor, used only be the makeVirtualSitesHandler() factory function
    VirtualSitesHandler(const gmx_mtop_t& mtop, gmx_domdec_t* domdec, PbcType pbcType);

    ~VirtualSitesHandler();

    //! Returns the number of virtual sites acting over multiple update groups
    int numInterUpdategroupVirtualSites() const;

    //! Set VSites and distribute VSite work over threads, should be called after each DD partitioning
    void setVirtualSites(ArrayRef<const InteractionList> ilist, const t_mdatoms& mdatoms);

    /*! \brief Create positions of vsite atoms based for the local system
     *
     * \param[in,out] x        The coordinates
     * \param[in]     dt       The time step
     * \param[in,out] v        When not empty, velocities for vsites are set as displacement/dt
     * \param[in]     box      The box
     */
    void construct(ArrayRef<RVec> x, real dt, ArrayRef<RVec> v, const matrix box) const;

    //! Tells how to handle virial contributions due to virtual sites
    enum class VirialHandling : int
    {
        None,     //!< Do not compute virial contributions
        Pbc,      //!< Add contributions working over PBC to shift forces
        NonLinear //!< Compute contributions due to non-linear virtual sites
    };

    /*! \brief Spread the force operating on the vsite atoms on the surrounding atoms.
     *
     * vsite should point to a valid object.
     * The virialHandling parameter determines how virial contributions are handled.
     * If this is set to Linear, shift forces are accumulated into fshift.
     * If this is set to NonLinear, non-linear contributions are added to virial.
     * This non-linear correction is required when the virial is not calculated
     * afterwards from the particle position and forces, but in a different way,
     * as for instance for the PME mesh contribution.
     */
    void spreadForces(ArrayRef<const RVec> x,
                      ArrayRef<RVec>       f,
                      VirialHandling       virialHandling,
                      ArrayRef<RVec>       fshift,
                      matrix               virial,
                      t_nrnb*              nrnb,
                      const matrix         box,
                      gmx_wallcycle*       wcycle);

private:
    //! Implementation type.
    class Impl;
    //! Implementation object.
    PrivateImplPointer<Impl> impl_;
};

/*! \brief Create positions of vsite atoms based for the local system
 *
 * \param[in,out] x        The coordinates
 * \param[in]     ip       Interaction parameters
 * \param[in]     ilist    The interaction list
 */
void constructVirtualSites(ArrayRef<RVec>                  x,
                           ArrayRef<const t_iparams>       ip,
                           ArrayRef<const InteractionList> ilist);

/*! \brief Create positions of vsite atoms for the whole system assuming all molecules are wholex
 *
 * \param[in]     mtop  The global topology
 * \param[in,out] x     The global coordinates
 */
void constructVirtualSitesGlobal(const gmx_mtop_t& mtop, ArrayRef<RVec> x);

//! Tells how to handle virial contributions due to virtual sites
enum class VirtualSiteVirialHandling : int
{
    None,     //!< Do not compute virial contributions
    Pbc,      //!< Add contributions working over PBC to shift forces
    NonLinear //!< Compute contributions due to non-linear virtual sites
};

//! Return the number of non-linear virtual site constructions in the system
int countNonlinearVsites(const gmx_mtop_t& mtop);

/*! \brief Return the number of virtual sites that cross update groups
 *
 * \param[in] mtop                           The global topology
 * \param[in] updateGroupingPerMoleculetype  Update grouping per molecule type, pass empty when not using update groups
 */
int countInterUpdategroupVsites(const gmx_mtop_t&                 mtop,
                                ArrayRef<const RangePartitioning> updateGroupingPerMoleculetype);

/*! \brief Create the virtual site handler
 *
 * \param[in] mtop      The global topology
 * \param[in] cr        The communication record
 * \param[in] pbcType   The type of PBC
 * \returns A valid vsite handler object or nullptr when there are no virtual sites
 */
std::unique_ptr<VirtualSitesHandler> makeVirtualSitesHandler(const gmx_mtop_t& mtop,
                                                             const t_commrec*  cr,
                                                             PbcType           pbcType);

} // namespace gmx

#endif
