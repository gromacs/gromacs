/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
//
// Created by Eric Irrgang on 11/10/17.
//

#include "gmxpre.h"

#include "restraintmdmodule.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"

#include "restraintmdmodule-impl.h"

gmx::RestraintForceProvider::RestraintForceProvider(std::shared_ptr<gmx::IRestraintPotential> restraint,
                                                    const std::vector<unsigned long int>     &sites) :
    restraint_ {std::move(restraint)},
sites_ {}
{
    assert(restraint_ != nullptr);
    assert(sites_.empty());
    for (auto && site : sites)
    {
        sites_.emplace_back(site);
    }
    assert(sites_.size() >= 2);
}

void gmx::RestraintForceProvider::calculateForces(const gmx::ForceProviderInput &forceProviderInput,
                                                  gmx::ForceProviderOutput     * forceProviderOutput)
{
    using gmx::detail::vec3;
    using gmx::detail::make_vec3;

    assert(restraint_ != nullptr);

    const auto &mdatoms {
        forceProviderInput.mdatoms_
    };
    assert(mdatoms.homenr >= 0);

    const auto &box {
        forceProviderInput.box_
    };
    assert(check_box(-1, box) == nullptr);
    t_pbc pbc {};
    set_pbc(&pbc, -1, box);

    const auto &x {
        forceProviderInput.x_
    };
    const auto &cr {
        forceProviderInput.cr_
    };
    const auto &t {
        forceProviderInput.t_
    };
    // Cooperatively get Cartesian coordinates for center of mass of each site
    RVec r1 = sites_[0].centerOfMass(cr,
                                     static_cast<size_t>(mdatoms.homenr),
                                     x,
                                     t);
    // r2 is to be constructed as
    // r2 = (site[N] - site[N-1]) + (site_{N-1} - site_{N-2}) + ... + (site_2 - site_1) + site_1
    // where the minimum image convention is applied to each path but not to the overall sum.
    // It is redundant to pass both r1 and r2 to called code, and potentially confusing, since
    // r1 may refer to an actual coordinate in the simulation while r2 may be in an expanded
    // Cartesian coordinate system. Called code should not use r1 and r2 to attempt to identify
    // sites in the simulation. If we need that functionality, we should do it separately by
    // allowing called code to look up atoms by tag or global index.
    RVec r2 {
        r1[0], r1[1], r1[2]
    };
    rvec dr {
        0, 0, 0
    };
    // Build r2 by following a path of difference vectors that are each presumed to be less than
    // a half-box apart, in case we are battling periodic boundary conditions along the lines of
    // a big molecule in a small box.
    for (size_t i = 0; i < sites_.size() - 1; ++i)
    {
        RVec a = sites_[i].centerOfMass(cr,
                                        static_cast<size_t>(mdatoms.homenr),
                                        x,
                                        t);
        RVec b = sites_[i + 1].centerOfMass(cr,
                                            static_cast<size_t>(mdatoms.homenr),
                                            x,
                                            t);
        // dr = minimum_image_vector(b - a)
        pbc_dx(&pbc, b, a, dr);
        r2[0] += dr[0];
        r2[1] += dr[1];
        r2[2] += dr[2];
    }
    // In the case of a single-atom site, r1 and r2 are now correct if local or [0,0,0] if not local.


    // Master rank update call-back. This needs to be moved to a discrete place in the
    // time step to avoid extraneous barriers. The code would be prettier with "futures"...
    if ((cr.dd == nullptr) || MASTER(&cr))
    {
        restraint_->update(make_vec3<real>(r1[0],
                                           r1[1],
                                           r1[2]),
                           make_vec3<real>(r2[0],
                                           r2[1],
                                           r2[2]),
                           t);
    }
    // All ranks wait for the update to finish.
    // tMPI ranks are depending on structures that may have just been updated.
    if (DOMAINDECOMP(&cr))
    {
        // Note: this assumes that all ranks are hitting this line, which is not generally true.
        // I need to find the right subcommunicator. What I really want is a _scoped_ communicator...
        gmx_barrier(&cr);
    }

    // Apply restraint on all thread ranks only after any updates have been made.
    auto result = restraint_->evaluate(make_vec3<real>(r1[0], r1[1], r1[2]),
                                       make_vec3<real>(r2[0], r2[1], r2[2]),
                                       t);

    // This can easily be generalized for pair restraints that apply to selections instead of
    // individual indices, or to restraints that aren't pair restraints.
    const int site1 {
        static_cast<int>(sites_.front().index())
    };
    const int* aLocal {
        &site1
    };
    // Set forces using index `site1` if no domain decomposition, otherwise set with local index if available.
    rvec *force = as_rvec_array(forceProviderOutput->forceWithVirial_.force_.data());
    if ((cr.dd == nullptr) || (aLocal = cr.dd->ga2la->findHome(site1)))
    {
        force[static_cast<size_t>(*aLocal)][0] += result.force.x;
        force[static_cast<size_t>(*aLocal)][1] += result.force.y;
        force[static_cast<size_t>(*aLocal)][2] += result.force.z;
    }

    // Note: Currently calculateForces is called once per restraint and each restraint
    // applies to a pair of atoms. Future optimizations may consolidate multiple restraints
    // with possibly duplicated sites, in which case we may prefer to iterate over non-frozen
    // sites to apply forces without explicitly expressing pairwise symmetry as in the
    // following logic.
    const int site2 {
        static_cast<int>(sites_.back().index())
    };
    const int* bLocal {
        &site2
    };
    if ((cr.dd == nullptr) || (bLocal = cr.dd->ga2la->findHome(site2)))
    {
        force[static_cast<size_t>(*bLocal)][0] -= result.force.x;
        force[static_cast<size_t>(*bLocal)][1] -= result.force.y;
        force[static_cast<size_t>(*bLocal)][2] -= result.force.z;
    }
}

gmx::RestraintMDModuleImpl::~RestraintMDModuleImpl() = default;

gmx::RestraintMDModuleImpl::RestraintMDModuleImpl(std::shared_ptr<gmx::IRestraintPotential> restraint,
                                                  const std::vector<unsigned long int>     &sites) :
    forceProvider_ {::gmx::compat::make_unique<RestraintForceProvider>(restraint, sites)},
outputProvider_ {
    ::gmx::compat::make_unique<RestraintOutputProvider>()
},
optionProvider_ {
    ::gmx::compat::make_unique<RestraintOptionProvider>()
}
{
    assert(forceProvider_ != nullptr);
    assert(outputProvider_ != nullptr);
    assert(optionProvider_ != nullptr);
}

gmx::IMdpOptionProvider *gmx::RestraintMDModuleImpl::mdpOptionProvider()
{
    assert(optionProvider_ != nullptr);
    return optionProvider_.get();
}

gmx::IMDOutputProvider *gmx::RestraintMDModuleImpl::outputProvider()
{
    assert(outputProvider_ != nullptr);
    return outputProvider_.get();
}

void gmx::RestraintMDModuleImpl::initForceProviders(ForceProviders *forceProviders)
{
    assert(forceProvider_ != nullptr);
    assert(forceProviders != nullptr);
    forceProviders->addForceProvider(forceProvider_.get());
}


// Needs to be defined after implementation type is complete in order to have unique_ptr member.
gmx::RestraintMDModule::~RestraintMDModule() = default;



gmx::IMdpOptionProvider *gmx::RestraintMDModule::mdpOptionProvider()
{
    assert(impl_ != nullptr);
    return impl_->mdpOptionProvider();
}

gmx::IMDOutputProvider *gmx::RestraintMDModule::outputProvider()
{
    assert(impl_ != nullptr);
    return impl_->outputProvider();
}

void gmx::RestraintMDModule::initForceProviders(ForceProviders *forceProviders)
{
    assert(impl_ != nullptr);
    impl_->initForceProviders(forceProviders);
}

std::unique_ptr<gmx::RestraintMDModule>
gmx::RestraintMDModule::create(std::shared_ptr<gmx::IRestraintPotential> restraint, const std::vector<unsigned long int> &sites)
{
    auto implementation = ::gmx::compat::make_unique<RestraintMDModuleImpl>(std::move(restraint), sites);
    std::unique_ptr<gmx::RestraintMDModule> newModule {
        new RestraintMDModule(std::move(implementation))
    };
    return newModule;
}

// private constructor to implement static create() method.
gmx::RestraintMDModule::RestraintMDModule(std::unique_ptr<RestraintMDModuleImpl> restraint) :
    impl_ {std::move(restraint)}
{
}
