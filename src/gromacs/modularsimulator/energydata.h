/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
 * \brief Declares the energy element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_ENERGYELEMENT_MICROSTATE_H
#define GMX_ENERGYELEMENT_MICROSTATE_H

#include "gromacs/mdtypes/state.h"

#include "modularsimulatorinterfaces.h"

class gmx_ekindata_t;
struct gmx_enerdata_t;
struct gmx_mtop_t;
struct ObservablesHistory;
struct pull_t;
struct t_fcdata;
struct t_inputrec;
struct SimulationGroups;

namespace gmx
{
enum class StartingBehavior;
class Constraints;
class EnergyOutput;
class FreeEnergyPerturbationData;
class GlobalCommunicationHelper;
class LegacySimulatorData;
class MDAtoms;
class ModularSimulatorAlgorithmBuilderHelper;
class ObservablesReducer;
class ParrinelloRahmanBarostat;
class StatePropagatorData;
class VelocityScalingTemperatureCoupling;
struct MDModulesNotifiers;

//! Function type for elements contributing energy
using EnergyContribution = std::function<real(Step, Time)>;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Data class managing energies
 *
 * The EnergyData owns the EnergyObject,
 * the tensors for the different virials and the pressure as well as
 * the total dipole vector. It has a member class which is part of the
 * simulator loop and and is responsible
 * for saving energy data and writing it to trajectory.
 *
 * The EnergyData offers an interface to add virial contributions,
 * but also allows access to the raw pointers to tensor data, the
 * dipole vector, and the legacy energy data structures.
 *
 * The EnergyData owns an object of type EnergyData::Element,
 * which takes part in the simulation loop, allowing to record
 * and output energies during the simulation.
 */
class EnergyData final
{
public:
    //! Constructor
    EnergyData(StatePropagatorData*        statePropagatorData,
               FreeEnergyPerturbationData* freeEnergyPerturbationData,
               const gmx_mtop_t&           globalTopology,
               const t_inputrec*           inputrec,
               const MDAtoms*              mdAtoms,
               gmx_enerdata_t*             enerd,
               gmx_ekindata_t*             ekind,
               const Constraints*          constr,
               FILE*                       fplog,
               t_fcdata*                   fcd,
               const MDModulesNotifiers&   mdModulesNotifiers,
               bool                        isMasterRank,
               ObservablesHistory*         observablesHistory,
               StartingBehavior            startingBehavior,
               bool                        simulationsShareState,
               pull_t*                     pullWork);

    /*! \brief Final output
     *
     * Prints the averages to log. This is called from ModularSimulatorAlgorithm.
     *
     * \see ModularSimulatorAlgorithm::teardown
     */
    void teardown();

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

    /*! \brief Get const pointer to energy structure
     *
     */
    const gmx_enerdata_t* enerdata() const;

    /*! \brief Get pointer to kinetic energy structure
     *
     */
    gmx_ekindata_t* ekindata();

    /*! \brief Get pointer to needToSumEkinhOld
     *
     */
    bool* needToSumEkinhOld();

    /*! \brief Whether kinetic energy was read from checkpoint
     *
     * This is needed by the compute globals element
     * TODO: Remove this when moving global reduction to client system (#3421)
     */
    [[nodiscard]] bool hasReadEkinFromCheckpoint() const;

    /*! \brief Add conserved energy contribution
     *
     * This allows other elements to register callbacks for contributions to
     * the conserved energy term.
     */
    void addConservedEnergyContribution(EnergyContribution&& energyContribution);

    /*! \brief set Parrinello-Rahman barostat
     *
     * This allows to set a pointer to the Parrinello-Rahman barostat used to
     * print the box velocities.
     */
    void setParrinelloRahmanBoxVelocities(std::function<const rvec*()>&& parrinelloRahmanBoxVelocities);

    /*! \brief Initialize energy history
     *
     * Kept as a static function to allow usage from legacy code
     * \todo Make member function once legacy use is not needed anymore
     */
    static void initializeEnergyHistory(StartingBehavior    startingBehavior,
                                        ObservablesHistory* observablesHistory,
                                        EnergyOutput*       energyOutput);

    /*! \brief Request (local) kinetic energy update
     */
    void updateKineticEnergy();

    //! The element taking part in the simulator loop
    class Element;
    //! Get pointer to element (whose lifetime is managed by this)
    Element* element();

private:
    /*! \brief Setup (needs file pointer)
     *
     * Initializes the EnergyOutput object, and does some logging output.
     *
     * \param mdoutf  File pointer
     */
    void setup(gmx_mdoutf* mdoutf);

    /*! \brief Save data at energy steps
     *
     * \param step  The current step
     * \param time  The current time
     * \param isEnergyCalculationStep  Whether the current step is an energy calculation step
     * \param isFreeEnergyCalculationStep  Whether the current step is a free energy calculation step
     */
    void doStep(Step step, Time time, bool isEnergyCalculationStep, bool isFreeEnergyCalculationStep);

    /*! \brief Write to energy trajectory
     *
     * This is only called by master - writes energy to trajectory and to log.
     */
    void write(gmx_mdoutf* outf, Step step, Time time, bool writeTrajectory, bool writeLog);

    /*
     * Data owned by EnergyData
     */
    //! The element
    std::unique_ptr<Element> element_;
    //! The energy output object
    std::unique_ptr<EnergyOutput> energyOutput_;
    //! Helper object to checkpoint kinetic energy data
    ekinstate_t ekinstate_;

    //! Whether this is the master rank
    const bool isMasterRank_;

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
    //! Whether we have read ekin from checkpoint
    bool hasReadEkinFromCheckpoint_;

    //! Describes how the simulation (re)starts
    const StartingBehavior startingBehavior_;

    /*
     * Pointers to Simulator data
     */
    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the state propagator data
    StatePropagatorData* statePropagatorData_;
    //! Pointer to the free energy perturbation data
    FreeEnergyPerturbationData* freeEnergyPerturbationData_;

    //! Callbacks contributing to the conserved energy term
    std::vector<EnergyContribution> conservedEnergyContributions_;
    //! Callback to the Parrinello-Rahman box velocities
    std::function<const rvec*()> parrinelloRahmanBoxVelocities_;

    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Full system topology.
    const gmx_mtop_t& top_global_;
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
    //! Notifiers to MD modules
    const MDModulesNotifiers& mdModulesNotifiers_;
    //! Global topology groups
    const SimulationGroups* groups_;
    //! History of simulation observables.
    ObservablesHistory* observablesHistory_;
    //! Whether simulations share the state
    bool simulationsShareState_;
    //! The pull work object.
    pull_t* pullWork_;
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Element for EnergyData
 *
 * This member class allows EnergyData to take part in the simulator
 * loop.
 *
 * It subscribes to the trajectory signaller, the energy signaller,
 * and the logging signaller to know when an energy calculation is
 * needed and when a non-recording step is enough. The simulator
 * builder is responsible to place the element in a location at
 * which a valid energy state is available. The EnergyData::Element is
 * also a subscriber to the trajectory writer element, as it is
 * responsible to write energy data to trajectory.
 */
class EnergyData::Element final :
    public ISimulatorElement,
    public ITrajectoryWriterClient,
    public ITrajectorySignallerClient,
    public IEnergySignallerClient,
    public ICheckpointHelperClient
{
public:
    //! Constructor
    Element(EnergyData* energyData, bool isMasterRank);

    /*! \brief Register run function for step / time
     *
     * This needs to be called when the energies are at a full time step.
     * Positioning this element is the responsibility of the programmer.
     *
     * This is also the place at which the current state becomes the previous
     * state.
     *
     * \param step                 The step number
     * \param time                 The time
     * \param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    //! No element setup needed
    void elementSetup() override {}

    //! No element teardown needed
    void elementTeardown() override {}

    //! ICheckpointHelperClient write checkpoint implementation
    void saveCheckpointState(std::optional<WriteCheckpointData> checkpointData, const t_commrec* cr) override;
    //! ICheckpointHelperClient read checkpoint implementation
    void restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData, const t_commrec* cr) override;
    //! ICheckpointHelperClient key implementation
    const std::string& clientID() override;

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
    EnergyData* energyData_;

    /*! \brief Setup (needs file pointer)
     *
     * ITrajectoryWriterClient implementation.
     *
     * Initializes the EnergyOutput object, and does some logging output.
     *
     * \param mdoutf  File pointer
     */
    void trajectoryWriterSetup(gmx_mdoutf* mdoutf) override;
    //! No trajectory writer teardown needed
    void trajectoryWriterTeardown(gmx_mdoutf gmx_unused* outf) override {}

    //! ITrajectoryWriterClient implementation.
    std::optional<SignallerCallback> registerTrajectorySignallerCallback(TrajectoryEvent event) override;
    //! ITrajectorySignallerClient implementation
    std::optional<ITrajectoryWriterCallback> registerTrajectoryWriterCallback(TrajectoryEvent event) override;
    //! IEnergySignallerClient implementation
    std::optional<SignallerCallback> registerEnergyCallback(EnergySignallerEvent event) override;


    //! CheckpointHelper identifier
    const std::string identifier_ = "EnergyElement";
    //! Helper function to read from / write to CheckpointData
    template<CheckpointDataOperation operation>
    void doCheckpointData(CheckpointData<operation>* checkpointData);

    //! Whether this is the master rank
    const bool isMasterRank_;
    //! The next communicated energy writing step
    Step energyWritingStep_;
    //! The next communicated energy calculation step
    Step energyCalculationStep_;
    //! The next communicated free energy calculation step
    Step freeEnergyCalculationStep_;
};

} // namespace gmx

#endif // GMX_ENERGYELEMENT_MICROSTATE_H
