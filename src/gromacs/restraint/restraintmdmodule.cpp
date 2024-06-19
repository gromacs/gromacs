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

#include "gmxpre.h"

#include "restraintmdmodule.h"

#include <cstddef>

#include <memory>
#include <utility>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "restraintmdmodule_impl.h"

namespace gmx
{
class IMDOutputProvider;
class IMdpOptionProvider;

RestraintForceProvider::RestraintForceProvider(std::shared_ptr<IRestraintPotential> restraint,
                                               const std::vector<int>&              sites) :
    restraint_{ std::move(restraint) }
{
    GMX_ASSERT(restraint_, "Valid RestraintForceProviders wrap non-null restraints.");
    GMX_ASSERT(sites_.empty(), "");
    for (auto&& site : sites)
    {
        sites_.emplace_back(site);
    }
    if (sites_.size() < 2)
    {
        GMX_THROW(InvalidInputError("Restraints require at least two sites to calculate forces."));
    }
}

void RestraintForceProvider::calculateForces(const ForceProviderInput& forceProviderInput,
                                             ForceProviderOutput*      forceProviderOutput)
{
    GMX_ASSERT(restraint_, "Restraint must be initialized.");

    const int homenr = forceProviderInput.homenr_;
    GMX_ASSERT(homenr >= 0, "number of home atoms must be non-negative.");

    const auto& box = forceProviderInput.box_;
    GMX_ASSERT(check_box(PbcType::Unset, box) == nullptr, "Invalid box.");
    t_pbc pbc{};
    set_pbc(&pbc, PbcType::Unset, box);

    const auto& x  = forceProviderInput.x_;
    const auto& cr = forceProviderInput.cr_;
    const auto& t  = forceProviderInput.t_;
    // Cooperatively get Cartesian coordinates for center of mass of each site
    RVec r1 = sites_[0].centerOfMass(cr, static_cast<size_t>(homenr), x, t);
    // r2 is to be constructed as
    // r2 = (site[N] - site[N-1]) + (site_{N-1} - site_{N-2}) + ... + (site_2 - site_1) + site_1
    // where the minimum image convention is applied to each path but not to the overall sum.
    // It is redundant to pass both r1 and r2 to called code, and potentially confusing, since
    // r1 may refer to an actual coordinate in the simulation while r2 may be in an expanded
    // Cartesian coordinate system. Called code should not use r1 and r2 to attempt to identify
    // sites in the simulation. If we need that functionality, we should do it separately by
    // allowing called code to look up atoms by tag or global index.
    RVec r2 = { r1[0], r1[1], r1[2] };
    rvec dr = { 0, 0, 0 };
    // Build r2 by following a path of difference vectors that are each presumed to be less than
    // a half-box apart, in case we are battling periodic boundary conditions along the lines of
    // a big molecule in a small box.
    for (size_t i = 0; i < sites_.size() - 1; ++i)
    {
        RVec a = sites_[i].centerOfMass(cr, static_cast<size_t>(homenr), x, t);
        RVec b = sites_[i + 1].centerOfMass(cr, static_cast<size_t>(homenr), x, t);
        // dr = minimum_image_vector(b - a)
        pbc_dx(&pbc, b, a, dr);
        r2[0] += dr[0];
        r2[1] += dr[1];
        r2[2] += dr[2];
    }
    // In the case of a single-atom site, r1 and r2 are now correct if local or [0,0,0] if not local.


    // Main rank update call-back. This needs to be moved to a discrete place in the
    // time step to avoid extraneous barriers. The code would be prettier with "futures"...
    if ((cr.dd == nullptr) || MAIN(&cr))
    {
        restraint_->update(RVec(r1), r2, t);
    }
    // All ranks wait for the update to finish.
    // tMPI ranks are depending on structures that may have just been updated.
    if (haveDDAtomOrdering(cr))
    {
        // Note: this assumes that all ranks are hitting this line, which is not generally true.
        // I need to find the right subcommunicator. What I really want is a _scoped_ communicator...
        gmx_barrier(cr.mpi_comm_mygroup);
    }

    // Apply restraint on all thread ranks only after any updates have been made.
    auto result = restraint_->evaluate(RVec(r1), r2, t);

    // This can easily be generalized for pair restraints that apply to selections instead of
    // individual indices, or to restraints that aren't pair restraints.
    const int  site1  = static_cast<int>(sites_.front().index());
    const int* aLocal = &site1;
    // Set forces using index `site1` if no domain decomposition, otherwise set with local index if available.
    const auto& force = forceProviderOutput->forceWithVirial_.force_;
    if ((cr.dd == nullptr) || (aLocal = cr.dd->ga2la->findHome(site1)))
    {
        force[static_cast<size_t>(*aLocal)] += result.force;
    }

    // Note: Currently calculateForces is called once per restraint and each restraint
    // applies to a pair of atoms. Future optimizations may consolidate multiple restraints
    // with possibly duplicated sites, in which case we may prefer to iterate over non-frozen
    // sites to apply forces without explicitly expressing pairwise symmetry as in the
    // following logic.
    const int  site2  = static_cast<int>(sites_.back().index());
    const int* bLocal = &site2;
    if ((cr.dd == nullptr) || (bLocal = cr.dd->ga2la->findHome(site2)))
    {
        force[static_cast<size_t>(*bLocal)] -= result.force;
    }
}

RestraintMDModuleImpl::~RestraintMDModuleImpl() = default;

RestraintMDModuleImpl::RestraintMDModuleImpl(std::shared_ptr<IRestraintPotential> restraint,
                                             const std::vector<int>&              sites) :
    forceProvider_(std::make_unique<RestraintForceProvider>(restraint, sites))
{
    GMX_ASSERT(forceProvider_, "Class invariant implies non-null ForceProvider.");
}

void RestraintMDModuleImpl::initForceProviders(ForceProviders* forceProviders)
{
    GMX_ASSERT(forceProvider_, "Class invariant implies non-null ForceProvider member.");
    GMX_ASSERT(forceProviders, "Provided ForceProviders* assumed to be non-null.");
    forceProviders->addForceProvider(forceProvider_.get());
}


// Needs to be defined after implementation type is complete in order to have unique_ptr member.
RestraintMDModule::~RestraintMDModule() = default;


IMdpOptionProvider* RestraintMDModule::mdpOptionProvider()
{
    return nullptr;
}

IMDOutputProvider* RestraintMDModule::outputProvider()
{
    return nullptr;
}

void RestraintMDModule::initForceProviders(ForceProviders* forceProviders)
{
    GMX_ASSERT(impl_, "Class invariant implies non-null implementation member.");
    impl_->initForceProviders(forceProviders);
}

std::unique_ptr<RestraintMDModule> RestraintMDModule::create(std::shared_ptr<IRestraintPotential> restraint,
                                                             const std::vector<int>& sites)
{
    auto implementation = std::make_unique<RestraintMDModuleImpl>(std::move(restraint), sites);
    auto newModule      = std::make_unique<RestraintMDModule>(std::move(implementation));
    return newModule;
}

void RestraintMDModule::subscribeToSimulationSetupNotifications(MDModulesNotifiers* /*notifiers*/)
{
}

void RestraintMDModule::subscribeToPreProcessingNotifications(MDModulesNotifiers* /*notifiers*/) {}

// private constructor to implement static create() method.
RestraintMDModule::RestraintMDModule(std::unique_ptr<RestraintMDModuleImpl> restraint) :
    impl_{ std::move(restraint) }
{
}

} // end namespace gmx
