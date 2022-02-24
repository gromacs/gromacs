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
 * \brief Declares the PME load balancing helper for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_PMELOADBALANCEHELPER_H
#define GMX_MODULARSIMULATOR_PMELOADBALANCEHELPER_H

#include "modularsimulatorinterfaces.h"

struct gmx_wallcycle;
struct pme_load_balancing_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;

namespace gmx
{
class MDLogger;
struct MdrunOptions;
class StatePropagatorData;
class SimulationWorkload;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Infrastructure element responsible for PME load balancing
 *
 * This encapsulates the function call to PME load balancing, which is
 * important for performance but outside of the current scope of the modular
 * simulator project. This relies on legacy data structures for the state.
 *
 * This element does not implement the ISimulatorElement interface, as
 * the Simulator is calling it explicitly between task queue population
 * steps. This allows elements to be aware of any changes before
 * deciding what functionality they need to run.
 */
class PmeLoadBalanceHelper final : public INeighborSearchSignallerClient
{
public:
    //! Constructor
    PmeLoadBalanceHelper(bool                 isVerbose,
                         StatePropagatorData* statePropagatorData,
                         FILE*                fplog,
                         t_commrec*           cr,
                         const MDLogger&      mdlog,
                         const t_inputrec*    inputrec,
                         gmx_wallcycle*       wcycle,
                         t_forcerec*          fr);

    //! Initialize the load balancing object
    void setup();
    //! Do load balancing
    void run(Step step, Time time);
    //! Final printout and deconstruction of the load balancing object
    void teardown();
    //! Whether PME load balancing printing is active \todo Check this!
    bool pmePrinting() const;

    //! Whether we're doing PME load balancing
    static bool doPmeLoadBalancing(const MdrunOptions&       mdrunOptions,
                                   const t_inputrec*         inputrec,
                                   const t_forcerec*         fr,
                                   const SimulationWorkload& simWorkload);

    //! Direct access to the load balancing object - used by reset counter
    const pme_load_balancing_t* loadBalancingObject();

private:
    //! The PME load balancing object - used by reset counter
    pme_load_balancing_t* pme_loadbal_;

    //! INeighborSearchSignallerClient implementation
    std::optional<SignallerCallback> registerNSCallback() override;

    //! The next NS step
    Step nextNSStep_;
    //! Whether we're being verbose
    const bool isVerbose_;
    //! Whether PME load balancing printing is active \todo Check this!
    bool bPMETunePrinting_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Handles communication.
    t_commrec* cr_;
    //! Handles logging.
    const MDLogger& mdlog_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wcycle_;
    //! Parameters for force calculations.
    t_forcerec* fr_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_PMELOADBALANCEHELPER_H
