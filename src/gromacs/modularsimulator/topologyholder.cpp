/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
#include "gromacs/modularsimulator/modularsimulatorinterfaces.h"
#include "gromacs/topology/topology.h"

namespace gmx
{
TopologyHolder::TopologyHolder(std::vector<ITopologyHolderClient*> clients,
                               const gmx_mtop_t&                   globalTopology,
                               gmx_localtop_t*                     localTopology,
                               const t_commrec*                    cr,
                               const t_inputrec*                   inputrec,
                               t_forcerec*                         fr,
                               MDAtoms*                            mdAtoms,
                               Constraints*                        constr,
                               VirtualSitesHandler*                vsite) :
    globalTopology_(globalTopology), localTopology_(localTopology), clients_(std::move(clients))
{
    if (!haveDDAtomOrdering(*cr))
    {
        // Generate and initialize new topology
        // Note that most of the data needed for the constructor is used here -
        // this function should probably be simplified sooner or later.
        // Note: Legacy mdrun resizes the force buffer in mdAlgorithmsSetupAtomData()
        //       TopologyHolder has no access to the forces, so we are passing a nullptr
        //       TODO: Find a unique approach to resizing the forces in modular simulator (#3461)
        mdAlgorithmsSetupAtomData(
                cr, *inputrec, globalTopology, localTopology_, fr, nullptr, mdAtoms, constr, vsite, nullptr);
    }
    // Send copy of initial topology to clients
    updateLocalTopology();
}

const gmx_mtop_t& TopologyHolder::globalTopology() const
{
    return globalTopology_;
}

void TopologyHolder::updateLocalTopology()
{
    for (auto& client : clients_)
    {
        client->setTopology(localTopology_);
    }
}
DomDecCallback TopologyHolder::registerDomDecCallback()
{
    return [this]() { updateLocalTopology(); };
}

void TopologyHolder::Builder::registerClient(ITopologyHolderClient* client)
{
    // Register client
    if (client)
    {
        if (state_ == ModularSimulatorBuilderState::NotAcceptingClientRegistrations)
        {
            throw SimulationAlgorithmSetupError(
                    "Tried to register to signaller after it was built.");
        }
        clients_.emplace_back(client);
    }
}
} // namespace gmx
