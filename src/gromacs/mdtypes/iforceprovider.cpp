/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements classes from iforceprovider.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "iforceprovider.h"

#include <utility>
#include <vector>

#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"

using namespace gmx;

class ForceProviders::Impl
{
public:
    Impl(gmx_wallcycle* wallCycle) : wallCycle_(wallCycle) {}

    std::vector<std::pair<IForceProvider*, std::optional<WallCycleCounter>>> providers_;
    gmx_wallcycle*                                                           wallCycle_;
};

ForceProviders::ForceProviders(gmx_wallcycle* wallCycle) : impl_(new Impl(wallCycle)) {}

ForceProviders::~ForceProviders() {}

void ForceProviders::addForceProvider(gmx::IForceProvider* provider, const std::string& cycleCounterName)
{
    std::optional<WallCycleCounter> counter;
    if (impl_->wallCycle_)
    {
        counter = impl_->wallCycle_->registerCycleCounter(cycleCounterName);
    }

    impl_->providers_.emplace_back(provider, counter);
}

bool ForceProviders::hasForceProvider() const
{
    return !impl_->providers_.empty();
}

void ForceProviders::calculateForces(const ForceProviderInput& forceProviderInput,
                                     ForceProviderOutput*      forceProviderOutput) const
{
    for (auto& provider : impl_->providers_)
    {
        if (provider.second)
        {
            wallcycle_start(impl_->wallCycle_, provider.second.value());
        }

        provider.first->calculateForces(forceProviderInput, forceProviderOutput);

        if (provider.second)
        {
            wallcycle_stop(impl_->wallCycle_, provider.second.value());
        }
    }
}
