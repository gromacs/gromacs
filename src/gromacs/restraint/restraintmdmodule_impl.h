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

#ifndef GROMACS_RESTRAINTMDMODULE_IMPL_H
#define GROMACS_RESTRAINTMDMODULE_IMPL_H

/*! \libinternal \file
 * \brief Implementation details for RestraintMDModule
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \ingroup module_restraint
 */

#include <iostream>
#include <mutex>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \libinternal
 * \brief Abstraction for a restraint interaction site.
 *
 * A restraint may operate on a single atom or some other entity, such as a selection of atoms.
 * The Restraint implementation is very independent from how coordinates are provided or what they mean.
 *
 * First implementation can only represent single atoms in a global context.
 *
 * Ultimately, this should be replaced with a more universal facility for acting
 * on distributed atom data or simple transformations thereof.
 *
 * \inlibraryapi
 * \ingroup module_restraint
 */
class Site
{
public:
    /*! \brief Construct from global atom indices
     *
     * \param globalIndex Atom index in the global state (as input to the simulation)
     */
    explicit Site(int globalIndex) : index_(globalIndex), r_(0, 0, 0) {}

    /*!
     * \brief Explicitly define copies.
     *
     * Implicit definition is not possible because of the mutex member,
     * and a copy constructor is necessary to use Site in a std::vector. Really, we should make
     * a copy point to the same implementation object to reuse its cache.

     */
    Site(const Site& site) : index_(site.index_), r_(site.r_) {}

    /*! \brief Disallow assignment.
     *
     * Assignment doesn't make sense because it implies that a site's meaning is fuzzy.
     * If the definition of a site is changing, just make a new site.
     * There's nothing to be gained by reusing one or by creating it uninitialized.
     */
    Site& operator=(const Site&) = delete;

    /*!
     * \brief Get the global atom index of an atomic site.
     *
     * \return global index provided at construction.
     *
     */
    int index() const { return index_; }

    /*!
     * \brief Get the position of this site at time t.
     *
     * \param cr Communications record.
     * \param nx Number of locally available atoms (size of local atom data arrays)
     * \param x Array of locally available atom coordinates.
     * \param t the current time.
     * \return position vector.
     *
     * \internal
     * By providing the current time, we can cache results in order to use them once per timestep.
     * In the long term, we would prefer to also allow client code to preregister interest in a
     * position at a given time, or issue "futures".
     */
    RVec centerOfMass(const t_commrec& cr, size_t nx, ArrayRef<const RVec> x, double gmx_unused t)
    {
        // Center of mass to return for the site. Currently the only form of site
        // implemented is as a global atomic coordinate.
        gmx::RVec r = { 0, 0, 0 };
        if (DOMAINDECOMP(&cr)) // Domain decomposition
        {
            // Get global-to-local indexing structure
            auto crossRef = cr.dd->ga2la;
            GMX_ASSERT(crossRef, "Domain decomposition must provide global/local cross-reference.");
            if (const auto localIndex = crossRef->findHome(index_))
            {
                GMX_ASSERT(localIndex,
                           "Expect not to reach this point if findHome does not find index_.");
                GMX_ASSERT(*localIndex < static_cast<decltype(*localIndex)>(nx),
                           "We assume that the local index cannot be larger than the number of "
                           "atoms.");
                GMX_ASSERT(*localIndex >= 0, "localIndex is a signed type, but is assumed >0.");
                // If atom is local, get its location
                copy_rvec(x[*localIndex], r);
            }
            else
            {
                // Nothing to contribute on this rank. Leave position == [0,0,0].
            }
            // AllReduce across the ranks of the simulation to get the center-of-mass
            // of the site locally available everywhere. For single-atom sites, this
            // is trivial: exactly one rank should have a non-zero position.
            // For future multi-atom selections,
            // we will receive weighted center-of-mass contributions from
            // each rank and combine to get the global center of mass.
            // \todo use generalized "pull group" facility when available.
            std::array<double, 3> buffer{ { r[0], r[1], r[2] } };
            // This should be an all-reduce sum, which gmx_sumd appears to be.
            gmx_sumd(3, buffer.data(), &cr);
            r[0] = static_cast<real>(buffer[0]);
            r[1] = static_cast<real>(buffer[1]);
            r[2] = static_cast<real>(buffer[2]);

        } // end domain decomposition branch
        else
        {
            // No DD so all atoms are local.
            copy_rvec(x[index_], r);
            (void)nx;
        }
        // Update cache and cache status.
        copy_rvec(r, r_);

        return r_;
    }

private:
    /*!
     * \brief Global index of the single-atom site.
     *
     * \todo This class should be a specialization of a more general Site data source.
     * \todo use LocalAtomSet
     */
    const int index_;

    /*!
     * \brief Last known value of the center-of-mass.
     *
     * Updated with centerOfMass().
     */
    RVec r_;
};

/*! \internal
 * \brief Provide IForceProvider for RestraintMDModuleImpl
 *
 * Adapter class from IForceProvider to IRestraintPotential.
 * Objects of this type are uniquely owned by instances of RestraintMDModuleImpl. The object will
 * dispatch calls to IForceProvider->calculateForces() to the functor managed by
 * RestraintMDModuleImpl. \ingroup module_restraint
 */
class RestraintForceProvider final : public gmx::IForceProvider
{
public:
    /*!
     * \brief Can only be constructed when initialized from a restraint.
     */
    RestraintForceProvider() = delete;

    ~RestraintForceProvider() = default;

    /*!
     * \brief RAII construction with an IRestraintPotential
     *
     * Note, this object must outlive the pointer that will be provided to ForceProviders.
     * \param restraint handle to an object providing restraint potential calculation
     * \param sites List of atomic site indices
     */
    explicit RestraintForceProvider(std::shared_ptr<gmx::IRestraintPotential> restraint,
                                    const std::vector<int>&                   sites);

    /*!
     * \brief Implement the IForceProvider interface.
     *
     * Update the force array with restraint contribution(s) for local atoms.
     *
     * RestraintForceProvider is implemented with the assumption that few
     * restraints apply to many atoms.
     * That is, the number of restraints affecting a large number of atoms is small,
     * though there may be several restraints that apply to few atoms each.
     * Under this assumption, it is considered computationally inexpensive to iterate
     * over restraints in an outer loop and iterate over atoms within each restraint.
     * This would be an invalid assumption if, say, several restraints applied
     * to an entire membrane or the entire solvent group.
     *
     * If the assumption causes performance problems, we can look for a good
     * way to reduce from several restraints in a single pass or a very
     * lightweight way to determine whether a given restraint applies to a given atom.
     * There is also the notion in the pulling code of a limited number of
     * "pull groups" used by the "pull coordinates".
     * The right optimization will depend on how the code is being used.
     *
     * Call the evaluator(s) for the restraints for the configured sites.
     * Forces are applied to atoms in the first and last site listed.
     * Intermediate sites are used as reference coordinates when the relevant
     * vector between sites is on the order of half a box length or otherwise
     * ambiguous in the case of periodic boundary conditions.
     */
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;

private:
    std::shared_ptr<gmx::IRestraintPotential> restraint_;
    std::vector<Site>                         sites_;
};

/*! \internal
 * \brief IMDModule implementation for RestraintMDModule.
 *
 * Provides IMDModule interface.
 *
 * \ingroup module_restraint
 */
class RestraintMDModuleImpl final
{
public:
    RestraintMDModuleImpl() = delete;
    /*!
     * \brief Wrap an object implementing IRestraintPotential
     *
     * \param restraint handle to restraint to wrap.
     * \param sites list of sites for framework to process for restraint force calculator.
     */
    RestraintMDModuleImpl(std::shared_ptr<gmx::IRestraintPotential> restraint,
                          const std::vector<int>&                   sites);

    /*!
     * \brief Allow moves.
     *
     * \{
     */
    RestraintMDModuleImpl(RestraintMDModuleImpl&&) noexcept = default;
    RestraintMDModuleImpl& operator=(RestraintMDModuleImpl&&) noexcept = default;
    /*! \} */

    ~RestraintMDModuleImpl();

    /*!
     * \brief Unused implementation of IMDModule aspect
     */
    IMdpOptionProvider* mdpOptionProvider();

    /*!
     * \brief Unused implementation of IMDModule aspect
     */
    IMDOutputProvider* outputProvider();

    /*!
     * \brief Implement IMDModule interface.
     *
     * \param forceProviders force module manager in the force record that will call this.
     *
     * The calling code must ensure that this object stays alive as long as forceProviders needs
     * the RestraintForceProvider, since forceProviders can't. Typically that is the duration of a do_md() call.
     */
    void initForceProviders(ForceProviders* forceProviders);

    //! handle to RestraintForceProvider implementation
    std::unique_ptr<RestraintForceProvider> forceProvider_;
};

} // end namespace gmx

#endif // GROMACS_RESTRAINTMDMODULE_IMPL_H
