/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief Declares the energy element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#ifndef GMX_ENERGYELEMENT_MICROSTATE_H
#define GMX_ENERGYELEMENT_MICROSTATE_H

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state.h"

#include "modularsimulatorinterfaces.h"

struct gmx_ekindata_t;
struct gmx_enerdata_t;
struct gmx_mtop_t;
struct ObservablesHistory;
struct t_fcdata;
struct t_inputrec;
struct SimulationGroups;

namespace gmx
{
enum class StartingBehavior;
class Constraints;
class EnergyOutput;
class FreeEnergyPerturbationElement;
class MDAtoms;
class ParrinelloRahmanBarostat;
class StatePropagatorData;
class VRescaleThermostat;
struct MdModulesNotifier;

/*! \libinternal
 * \ingroup module_modularsimulator
 * \brief Element managing energies
 *
 * The EnergyElement owns the EnergyObject, and is hence responsible
 * for saving energy data and writing it to trajectory. It also owns
 * the tensors for the different virials and the pressure as well as
 * the total dipole vector.
 *
 * It subscribes to the trajectory signaller, the energy signaller,
 * and the logging signaller to know when an energy calculation is
 * needed and when a non-recording step is enough. The simulator
 * builder is responsible to place the element in a location at
 * which a valid energy state is available. The EnergyElement is
 * also a subscriber to the trajectory writer element, as it is
 * responsible to write energy data to trajectory.
 *
 * The EnergyElement offers an interface to add virial contributions,
 * but also allows access to the raw pointers to tensor data, the
 * dipole vector, and the legacy energy data structures.
 */
class EnergyElement final :
    public ISimulatorElement,
    public ITrajectoryWriterClient,
    public ITrajectorySignallerClient,
    public IEnergySignallerClient,
    public ICheckpointHelperClient
{
public:
    //! Constructor
    EnergyElement(StatePropagatorData*           statePropagatorData,
                  FreeEnergyPerturbationElement* freeEnergyPerturbationElement,
                  const gmx_mtop_t*              globalTopology,
                  const t_inputrec*              inputrec,
                  const MDAtoms*                 mdAtoms,
                  gmx_enerdata_t*                enerd,
                  gmx_ekindata_t*                ekind,
                  const Constraints*             constr,
                  FILE*                          fplog,
                  t_fcdata*                      fcd,
                  const MdModulesNotifier&       mdModulesNotifier,
                  bool                           isMasterRank,
                  ObservablesHistory*            observablesHistory,
                  StartingBehavior               startingBehavior);

    /*! \brief Register run function for step / time
     *
     * This needs to be called when the energies are at a full time step.
     * Positioning this element is the responsibility of the programmer.
     *
     * This is also the place at which the current state becomes the previous
     * state.
     *
     * @param step                 The step number
     * @param time                 The time
     * @param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunctionPtr& registerRunFunction) override;

    //! No element setup needed
    void elementSetup() override {}

    /*! \brief Final output
     *
     * Prints the averages to log.
     */
    void elementTeardown() override;

    /*! \brief Add contribution to force virial
     *
     * This automatically resets the tensor if the step is higher
     * than the current step, starting the tensor calculation for
     * a new step at zero. Otherwise, it adds the new contribution
     * to the existing virial.
     */
    void addToForceVirial(const tensor virial, Step step);

    /*! \brief Add contribution to constraint virial
     *
     * This automatically resets the tensor if the step is higher
     * than the current step, starting the tensor calculation for
     * a new step at zero. Otherwise, it adds the new contribution
     * to the existing virial.
     */
    void addToConstraintVirial(const tensor virial, Step step);

    /*! \brief Get pointer to force virial tensor
     *
     * Allows access to the raw pointer to the tensor.
     */
    rvec* forceVirial(Step step);

    /*! \brief Get pointer to constraint virial tensor
     *
     * Allows access to the raw pointer to the tensor.
     */
    rvec* constraintVirial(Step step);

    /*! \brief Get pointer to total virial tensor
     *
     * Allows access to the raw pointer to the tensor.
     */
    rvec* totalVirial(Step step);

    /*! \brief Get pointer to pressure tensor
     *
     * Allows access to the raw pointer to the tensor.
     */
    rvec* pressure(Step step);

    /*! \brief Get pointer to mu_tot
     *
     * Allows access to the raw pointer to the dipole vector.
     */
    real* muTot();

    /*! \brief Get pointer to energy structure
     *
     */
    gmx_enerdata_t* enerdata();

    /*! \brief Get pointer to kinetic energy structure
     *
     */
    gmx_ekindata_t* ekindata();

    /*! \brief Get pointer to needToSumEkinhOld
     *
     */
    bool* needToSumEkinhOld();

    /*! \brief set vrescale thermostat
     *
     * This allows to set a pointer to the vrescale thermostat used to
     * print the thermostat integral.
     * TODO: This should be made obsolete my a more modular energy element
     */
    void setVRescaleThermostat(const VRescaleThermostat* vRescaleThermostat);

    /*! \brief set Parrinello-Rahman barostat
     *
     * This allows to set a pointer to the Parrinello-Rahman barostat used to
     * print the box velocities.
     * TODO: This should be made obsolete my a more modular energy element
     */
    void setParrinelloRahamnBarostat(const ParrinelloRahmanBarostat* parrinelloRahmanBarostat);

    /*! \brief Initialize energy history
     *
     * Kept as a static function to allow usage from legacy code
     * \todo Make member function once legacy use is not needed anymore
     */
    static void initializeEnergyHistory(StartingBehavior    startingBehavior,
                                        ObservablesHistory* observablesHistory,
                                        EnergyOutput*       energyOutput);

private:
    /*! \brief Setup (needs file pointer)
     *
     * ITrajectoryWriterClient implementation.
     *
     * Initializes the EnergyOutput object, and does some logging output.
     *
     * @param mdoutf  File pointer
     */
    void trajectoryWriterSetup(gmx_mdoutf* mdoutf) override;
    //! No trajectory writer teardown needed
    void trajectoryWriterTeardown(gmx_mdoutf gmx_unused* outf) override {}

    //! ITrajectoryWriterClient implementation.
    SignallerCallbackPtr registerTrajectorySignallerCallback(TrajectoryEvent event) override;
    //! ITrajectorySignallerClient implementation
    ITrajectoryWriterCallbackPtr registerTrajectoryWriterCallback(TrajectoryEvent event) override;
    //! IEnergySignallerClient implementation
    SignallerCallbackPtr registerEnergyCallback(EnergySignallerEvent event) override;

    /*! \brief Save data at energy steps
     *
     * @param time  The current time
     * @param isEnergyCalculationStep  Whether the current step is an energy calculation step
     * @param isFreeEnergyCalculationStep  Whether the current step is a free energy calculation step
     */
    void doStep(Time time, bool isEnergyCalculationStep, bool isFreeEnergyCalculationStep);

    /*! \brief Write to energy trajectory
     *
     * This is only called by master - writes energy to trajectory and to log.
     */
    void write(gmx_mdoutf* outf, Step step, Time time, bool writeTrajectory, bool writeLog);

    //! ICheckpointHelperClient implementation
    void writeCheckpoint(t_state* localState, t_state* globalState) override;

    /*
     * Data owned by EnergyElement
     */
    //! The energy output object
    std::unique_ptr<EnergyOutput> energyOutput_;

    //! Whether this is the master rank
    const bool isMasterRank_;
    //! The next communicated energy writing step
    Step energyWritingStep_;
    //! The next communicated energy calculation step
    Step energyCalculationStep_;
    //! The next communicated free energy calculation step
    Step freeEnergyCalculationStep_;

    //! The force virial tensor
    tensor forceVirial_;
    //! The constraint virial tensor
    tensor shakeVirial_;
    //! The total virial tensor
    tensor totalVirial_;
    //! The pressure tensor
    tensor pressure_;
    //! The total dipole moment
    rvec muTot_;

    //! The step number of the current force virial tensor
    Step forceVirialStep_;
    //! The step number of the current constraint virial tensor
    Step shakeVirialStep_;
    //! The step number of the current total virial tensor
    Step totalVirialStep_;
    //! The step number of the current pressure tensor
    Step pressureStep_;

    //! Whether ekinh_old needs to be summed up (set by compute globals)
    bool needToSumEkinhOld_;

    //! Describes how the simulation (re)starts
    const StartingBehavior startingBehavior_;

    //! Legacy state object used to communicate with energy output
    t_state dummyLegacyState_;

    /*
     * Pointers to Simulator data
     */
    //! Pointer to the state propagator data
    StatePropagatorData* statePropagatorData_;
    //! Pointer to the free energy perturbation element
    FreeEnergyPerturbationElement* freeEnergyPerturbationElement_;
    //! Pointer to the vrescale thermostat
    const VRescaleThermostat* vRescaleThermostat_;
    //! Pointer to the Parrinello-Rahman barostat
    const ParrinelloRahmanBarostat* parrinelloRahmanBarostat_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Full system topology.
    const gmx_mtop_t* top_global_;
    //! Atom parameters for this domain.
    const MDAtoms* mdAtoms_;
    //! Energy data structure
    gmx_enerdata_t* enerd_;
    //! Kinetic energy data
    gmx_ekindata_t* ekind_;
    //! Handles constraints.
    const Constraints* constr_;
    //! Handles logging.
    FILE* fplog_;
    //! Helper struct for force calculations.
    t_fcdata* fcd_;
    //! Notification to MD modules
    const MdModulesNotifier& mdModulesNotifier_;
    //! Global topology groups
    const SimulationGroups* groups_;
    //! History of simulation observables.
    ObservablesHistory* observablesHistory_;
};

} // namespace gmx

#endif // GMX_ENERGYELEMENT_MICROSTATE_H
