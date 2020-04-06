/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
 * \brief Defines the topology class for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "topologyholder.h"

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/mdsetup.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/topology/topology.h"

namespace gmx
{
TopologyHolder::TopologyHolder(const gmx_mtop_t& globalTopology,
                               const t_commrec*  cr,
                               const t_inputrec* inputrec,
                               t_forcerec*       fr,
                               MDAtoms*          mdAtoms,
                               Constraints*      constr,
                               gmx_vsite_t*      vsite) :
    globalTopology_(globalTopology),
    localTopology_(std::make_unique<gmx_localtop_t>(globalTopology.ffparams))
{
    if (!DOMAINDECOMP(cr))
    {
        // Generate and initialize new topology
        // Note that most of the data needed for the constructor is used here -
        // this function should probably be simplified sooner or later.
        // Note: Legacy mdrun resizes the force buffer in mdAlgorithmsSetupAtomData()
        //       TopologyHolder has no access to the forces, so we are passing a nullptr
        //       TODO: Find a unique approach to resizing the forces in modular simulator (#3461)
        mdAlgorithmsSetupAtomData(cr, inputrec, globalTopology, localTopology_.get(), fr, nullptr,
                                  mdAtoms, constr, vsite, nullptr);
    }
}

const gmx_mtop_t& TopologyHolder::globalTopology() const
{
    return globalTopology_;
}

void TopologyHolder::updateLocalTopology()
{
    for (auto& client : clients_)
    {
        client->setTopology(localTopology_.get());
    }
}

void TopologyHolder::registerClient(ITopologyHolderClient* client)
{
    // Register client
    clients_.emplace_back(client);
    // Send copy of current topology
    client->setTopology(localTopology_.get());
}
} // namespace gmx
