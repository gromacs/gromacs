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
 * \brief Declares the topology class for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */


#ifndef GMX_MODULARSIMULATOR_TOPOLOGYHOLDER_H
#define GMX_MODULARSIMULATOR_TOPOLOGYHOLDER_H

#include <memory>
#include <utility>
#include <vector>

#include "gromacs/compat/pointers.h"

#include "modularsimulatorinterfaces.h"

struct gmx_localtop_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;

namespace gmx
{
class Constraints;
class MDAtoms;
class VirtualSitesHandler;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Object holding the topology
 *
 * Clients can register to get an updated local topology whenever there
 * is a change (infrequent, only due to domdec currently).
 */
class TopologyHolder final : public IDomDecHelperClient
{
public:
    //! Constructor
    TopologyHolder(std::vector<ITopologyHolderClient*> clients,
                   const gmx_mtop_t&                   globalTopology,
                   gmx_localtop_t*                     localTopology,
                   const t_commrec*                    cr,
                   const t_inputrec*                   inputrec,
                   t_forcerec*                         fr,
                   MDAtoms*                            mdAtoms,
                   Constraints*                        constr,
                   VirtualSitesHandler*                vsite);

    //! Get global topology
    const gmx_mtop_t& globalTopology() const;

    //! Callback on domain decomposition repartitioning
    DomDecCallback registerDomDecCallback() override;

    //! Allow domdec to access local topology directly
    friend class DomDecHelper;

    //! The builder
    class Builder;

private:
    //! Constant reference to the global topology
    const gmx_mtop_t& globalTopology_;
    //! Pointer to the currently valid local topology
    gmx_localtop_t* localTopology_;

    //! List of clients to be updated if local topology changes
    std::vector<ITopologyHolderClient*> clients_;

    //! Update local topology
    void updateLocalTopology();
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Builder for the topology holder
 */
class TopologyHolder::Builder
{
public:
    //! Register topology client
    void registerClient(ITopologyHolderClient* client);

    //! Build TopologyHolder
    template<typename... Args>
    std::unique_ptr<TopologyHolder> build(Args&&... args);

private:
    //! List of clients to be updated if local topology changes
    std::vector<ITopologyHolderClient*> clients_;
    //! The state of the builder
    ModularSimulatorBuilderState state_ = ModularSimulatorBuilderState::AcceptingClientRegistrations;
};

template<typename... Args>
std::unique_ptr<TopologyHolder> TopologyHolder::Builder::build(Args&&... args)
{
    state_ = ModularSimulatorBuilderState::NotAcceptingClientRegistrations;
    return std::make_unique<TopologyHolder>(std::move(clients_), std::forward<Args>(args)...);
}

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_TOPOLOGYHOLDER_H
