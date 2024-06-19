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
 * \brief Declares the global reduction element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_COMPUTEGLOBALSELEMENT_H
#define GMX_MODULARSIMULATOR_COMPUTEGLOBALSELEMENT_H

#include <cstdio>

#include <functional>
#include <memory>
#include <optional>

#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/vcm.h"

#include "energydata.h"
#include "modularsimulatorinterfaces.h"
#include "statepropagatordata.h"
#include "topologyholder.h"

struct gmx_global_stat;
struct gmx_wallcycle;
struct t_nrnb;
struct gmx_mtop_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;

namespace gmx
{
class FreeEnergyPerturbationData;
class LegacySimulatorData;
class MDAtoms;
class MDLogger;
class ObservablesReducer;
class Constraints;
class GlobalCommunicationHelper;
class ModularSimulatorAlgorithmBuilderHelper;

//! \addtogroup module_modularsimulator
//! \{

//! The different global reduction schemes we know about
enum class ComputeGlobalsAlgorithm
{
    LeapFrog,
    VelocityVerlet
};

//! The function type allowing to request a check of the number of bonded interactions
typedef std::function<void()> CheckBondedInteractionsCallback;

/*! \internal
 * \brief Encapsulate the calls to `compute_globals`
 *
 * This element aims at offering an interface to the legacy
 * implementation which is compatible with the new simulator approach.
 *
 * The element comes in 3 (templated) flavors: the leap-frog case, the first
 * call during a velocity-verlet integrator, and the second call during a
 * velocity-verlet integrator. In velocity verlet, the state at the beginning
 * of the step corresponds to
 *     positions at time t
 *     velocities at time t - dt/2
 * The first velocity propagation (+dt/2) therefore actually corresponds to the
 * previous step, bringing the state to the full timestep at time t. Most global
 * reductions are made at this point. The second call is needed to correct the
 * constraint virial after the second propagation of velocities (+dt/2) and of
 * the positions (+dt).
 *
 * \tparam algorithm  The global reduction scheme
 */
template<ComputeGlobalsAlgorithm algorithm>
class ComputeGlobalsElement final : public ISimulatorElement, public IEnergySignallerClient, public ITrajectorySignallerClient
{
public:
    //! Constructor
    ComputeGlobalsElement(StatePropagatorData*        statePropagatorData,
                          EnergyData*                 energyData,
                          FreeEnergyPerturbationData* freeEnergyPerturbationData,
                          SimulationSignals*          signals,
                          int                         nstglobalcomm,
                          FILE*                       fplog,
                          const MDLogger&             mdlog,
                          t_commrec*                  cr,
                          const t_inputrec*           inputrec,
                          const MDAtoms*              mdAtoms,
                          t_nrnb*                     nrnb,
                          gmx_wallcycle*              wcycle,
                          t_forcerec*                 fr,
                          const gmx_mtop_t&           global_top,
                          Constraints*                constr,
                          ObservablesReducer*         observablesReducer);

    //! Destructor
    ~ComputeGlobalsElement() override;

    /*! \brief Element setup - first call to compute_globals
     *
     */
    void elementSetup() override;

    /*! \brief Register run function for step / time
     *
     * This registers the call to compute_globals when needed.
     *
     * \param step                 The step number
     * \param time                 The time
     * \param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    //! No element teardown needed
    void elementTeardown() override {}

    /*! \brief Factory method implementation
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper   Pointer to the \c GlobalCommunicationHelper object
     * \param observablesReducer          Pointer to the \c ObservablesReducer object
     *
     * \throws std::bad_any_cast  on internal error in VelocityVerlet algorithm builder.
     * \throws std::bad_alloc  when out of memory.
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement* getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                    StatePropagatorData*        statePropagatorData,
                                                    EnergyData*                 energyData,
                                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                    GlobalCommunicationHelper* globalCommunicationHelper,
                                                    ObservablesReducer*        observablesReducer);

private:
    //! IEnergySignallerClient implementation
    std::optional<SignallerCallback> registerEnergyCallback(EnergySignallerEvent event) override;
    //! ITrajectorySignallerClient implementation
    std::optional<SignallerCallback> registerTrajectorySignallerCallback(TrajectoryEvent event) override;
    //! The compute_globals call
    void compute(Step step, unsigned int flags, SimulationSignaller* signaller, bool useLastBox, bool isInit = false);

    //! Next step at which energy needs to be reduced
    Step energyReductionStep_;
    //! Next step at which virial needs to be reduced
    Step virialReductionStep_;

    //! For VV only, we need to schedule twice per step. This keeps track of the scheduling stage.
    Step vvSchedulingStep_;

    //! Whether center of mass motion stopping is enabled
    const bool doStopCM_;
    //! Number of steps after which center of mass motion is removed
    int nstcomm_;
    //! Compute globals communication period
    int nstglobalcomm_;
    //! The last (planned) step (only used for LF)
    const Step lastStep_;
    //! The initial step (only used for VV)
    const Step initStep_;
    //! A dummy signaller (used for setup and VV)
    std::unique_ptr<SimulationSignaller> nullSignaller_;

    //! Global reduction struct
    gmx_global_stat* gstat_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the microstate
    StatePropagatorData* statePropagatorData_;
    //! Pointer to the energy data (needed for the tensors and mu_tot)
    EnergyData* energyData_;
    //! Pointer to the free energy perturbation data
    FreeEnergyPerturbationData* freeEnergyPerturbationData_;

    //! Center of mass motion removal
    t_vcm vcm_;
    //! Signals
    SimulationSignals* signals_;

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Handles logging.
    const MDLogger& mdlog_;
    //! Handles communication.
    t_commrec* cr_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Full system topology - only needed for checkNumberOfBondedInteractions.
    const gmx_mtop_t& top_global_;
    //! Atom parameters for this domain.
    const MDAtoms* mdAtoms_;
    //! Handles constraints.
    Constraints* constr_;
    //! Manages flop accounting.
    t_nrnb* nrnb_;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wcycle_;
    //! Parameters for force calculations.
    t_forcerec* fr_;
    //! Coordinates reduction for observables
    ObservablesReducer* observablesReducer_;
};

//! \}
} // namespace gmx

#endif // GMX_MODULARSIMULATOR_COMPUTEGLOBALSELEMENT_H
