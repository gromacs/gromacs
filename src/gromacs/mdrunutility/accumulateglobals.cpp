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
/*! \internal \file
 * \brief Defines classes for components that collaborate to
 * accumulate simulation variables across PP ranks during an MD step.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include "accumulateglobals.h"

#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

constexpr bool c_debugBuild =
#ifndef NDEBUG
    true;
#else
    false;
#endif

void
AccumulateGlobalsBuilder::registerClient(IAccumulateGlobalsClient *newClient)
{
    if (newClient == nullptr)
    {
        return;
    }
    clients_.push_back(newClient);
}

AccumulateGlobals
AccumulateGlobalsBuilder::build() const
{
    AccumulateGlobals globals;

    // Ask each registered client how many doubles they require to be
    // reduced.
    std::vector<int>  globalsForEachClient;
    globalsForEachClient.reserve(clients_.size());
    size_t            numGlobals = 0;
    for (const auto &client : clients_)
    {
        globalsForEachClient.push_back(client->getNumGlobalsRequired());
        numGlobals += globalsForEachClient.back();
    }

    // Build the container for values to be reduced
    globals.globals_.resize(numGlobals);

    // Let the clients know where to find their globals.
    ArrayRef<double> globalsToSend = globals.globals_;
    size_t           currentStart  = 0;
    for (size_t i = 0; i != clients_.size(); ++i)
    {
        size_t numGlobalsForThisClient = globalsForEachClient[i];
        clients_[i]->setViewForGlobals(&globals, globalsToSend.subArray(currentStart, numGlobalsForThisClient));
        currentStart += numGlobalsForThisClient;
    }
    GMX_RELEASE_ASSERT(currentStart == numGlobals, "Failed to build current AccumulateGlobals");

    // Give the container a way to notify the clients that
    // there is data now, e.g. in debug mode.
    globals.clients_ = makeConstArrayRef(clients_);

    return globals;
}

void AccumulateGlobals::notifyReductionRequired()
{
    if (!reductionRequired_)
    {
        reductionRequired_ = true;
    }
}

bool AccumulateGlobals::isReductionRequired() const
{
    return reductionRequired_;
}

ArrayRef<double> AccumulateGlobals::getReductionView()
{
    return globals_;
}

void AccumulateGlobals::notifyClientsAfterCommunication()
{
    if (c_debugBuild)
    {
        for (const auto &client : clients_)
        {
            client->notifyAfterCommunication();
        }
    }
    reductionRequired_ = false;
}

} // namespace
