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
 * \brief Declares the domain decomposition helper for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_DOMDECHELPER_H
#define GMX_MODULARSIMULATOR_DOMDECHELPER_H

#include "modularsimulatorinterfaces.h"

struct gmx_localtop_t;
struct gmx_wallcycle;
struct pull_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_nrnb;

namespace gmx
{
class Constraints;
class FreeEnergyPerturbationData;
class ImdSession;
class MDAtoms;
class MDLogger;
class StatePropagatorData;
class TopologyHolder;
class VirtualSitesHandler;

//! \addtogroup module_modularsimulator
//! \{

/*! \internal
 * \brief Infrastructure element responsible for domain decomposition
 *
 * This encapsulates the function call to domain decomposition, which is
 * important for performance but outside of the current scope of the modular
 * simulator project. This relies on legacy data structures for the state
 * and the topology.
 *
 * This element does not implement the ISimulatorElement interface, as
 * the Simulator is calling it explicitly between task queue population
 * steps. This allows elements to be aware of any changes before
 * deciding what functionality they need to run.
 */
class DomDecHelper final : public INeighborSearchSignallerClient
{
public:
    //! Constructor
    DomDecHelper(bool                          isVerbose,
                 int                           verbosePrintInterval,
                 StatePropagatorData*          statePropagatorData,
                 TopologyHolder*               topologyHolder,
                 int                           nstglobalcomm,
                 FILE*                         fplog,
                 t_commrec*                    cr,
                 const MDLogger&               mdlog,
                 Constraints*                  constr,
                 const t_inputrec*             inputrec,
                 MDAtoms*                      mdAtoms,
                 t_nrnb*                       nrnb,
                 gmx_wallcycle*                wcycle,
                 t_forcerec*                   fr,
                 VirtualSitesHandler*          vsite,
                 ImdSession*                   imdSession,
                 pull_t*                       pull_work,
                 std::vector<DomDecCallback>&& domdecCallbacks);

    /*! \brief Run domain decomposition
     *
     * Does domain decomposition partitioning at neighbor searching steps
     *
     * \param step  The step number
     * \param time  The time
     */
    void run(Step step, Time time);

    /*! \brief The initial domain decomposition partitioning
     *
     */
    void setup();

private:
    //! INeighborSearchSignallerClient implementation
    std::optional<SignallerCallback> registerNSCallback() override;

    //! The next NS step
    Step nextNSStep_;
    //! Whether we're being verbose
    const bool isVerbose_;
    //! If we're being verbose, how often?
    const int verbosePrintInterval_;
    //! The global communication frequency
    const int nstglobalcomm_;

    //! List of callbacks to inform clients that DD happened
    std::vector<DomDecCallback> domdecCallbacks_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;
    //! Pointer to the topology
    TopologyHolder* topologyHolder_;

    //! Helper function unifying the DD partitioning calls in setup() and run()
    void partitionSystem(bool           verbose,
                         bool           isMasterState,
                         int            nstglobalcomm,
                         gmx_wallcycle* wcycle,
                         t_state*       localState,
                         t_state*       globalState);

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Handles communication.
    t_commrec* cr_;
    //! Handles logging.
    const MDLogger& mdlog_;
    //! Handles constraints.
    Constraints* constr_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Atom parameters for this domain.
    MDAtoms* mdAtoms_;
    //! Manages flop accounting.
    t_nrnb* nrnb_;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wcycle_;
    //! Parameters for force calculations.
    t_forcerec* fr_;
    //! Handles virtual sites.
    VirtualSitesHandler* vsite_;
    //! The Interactive Molecular Dynamics session.
    ImdSession* imdSession_;
    //! The pull work object.
    pull_t* pull_work_;
};

/*! \internal
 * \brief Builder for DomDecHelper
 *
 * This builder allows clients to register to the DomDecHelper in order to get
 * informed whenever system re-partitioning is performed.
 */
class DomDecHelperBuilder
{
public:
    //! Register DomDecHelper client
    void registerClient(IDomDecHelperClient* client);

    //! Return DomDecHelper instance
    template<typename... Args>
    std::unique_ptr<DomDecHelper> build(Args&&... args)
    {
        state_ = ModularSimulatorBuilderState::NotAcceptingClientRegistrations;
        std::vector<DomDecCallback> callbacks;
        for (const auto& client : clients_)
        {
            callbacks.emplace_back(client->registerDomDecCallback());
        }
        return std::make_unique<DomDecHelper>(std::forward<Args>(args)..., std::move(callbacks));
    }

private:
    //! List of clients to be updated after system partition
    std::vector<IDomDecHelperClient*> clients_;
    //! The state of the builder
    ModularSimulatorBuilderState state_ = ModularSimulatorBuilderState::AcceptingClientRegistrations;
};

//! \}
} // namespace gmx

#endif // GMX_MODULARSIMULATOR_DOMDECHELPER_H
