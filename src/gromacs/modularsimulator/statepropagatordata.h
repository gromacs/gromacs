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
 * \brief Declares the state for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_STATEPROPAGATORDATA_H
#define GMX_MODULARSIMULATOR_STATEPROPAGATORDATA_H

#include <cstdio>

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/real.h"

#include "modularsimulatorinterfaces.h"
#include "topologyholder.h"

struct gmx_mdoutf;
enum class PbcType : int;
struct t_commrec;
struct t_inputrec;
class t_state;
struct t_mdatoms;
struct t_trxframe;
struct gmx_mtop_t;

namespace gmx
{
enum class CheckpointDataOperation;
enum class ConstraintVariable;
class EnergyData;
class FreeEnergyPerturbationData;
class GlobalCommunicationHelper;
class LegacySimulatorData;
class ModularSimulatorAlgorithmBuilderHelper;
class ObservablesReducer;
enum class ReferenceTemperatureChangeAlgorithm;
template<CheckpointDataOperation operation>
class CheckpointData;
template<typename T>
class ArrayRef;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief StatePropagatorData and associated data
 *
 * The `StatePropagatorData` contains a little more than the pure
 * statistical-physical micro state, namely the positions,
 * velocities, forces, and box matrix, as well as a backup of
 * the positions and box of the last time step. While it takes
 * part in the simulator loop via its member class `Element`
 * to be able to backup positions /
 * boxes and save the current state if needed, it's main purpose
 * is to offer access to its data via getter methods. All elements
 * reading or writing to this data need a pointer to the
 * `StatePropagatorData` and need to request their data explicitly. This
 * will later simplify the understanding of data dependencies
 * between elements.
 *
 * Note that the `StatePropagatorData` can be converted to and from the
 * legacy `t_state` object. This is useful when dealing with
 * functionality which has not yet been adapted to use the new
 * data approach - of the elements currently implemented, only
 * domain decomposition, PME load balancing, and the initial
 * constraining are using this.
 */
class StatePropagatorData final
{
public:
    //! Constructor
    StatePropagatorData(int                numAtoms,
                        FILE*              fplog,
                        const t_commrec*   cr,
                        t_state*           globalState,
                        t_state*           localState,
                        bool               useGPU,
                        bool               canMoleculesBeDistributedOverPBC,
                        bool               writeFinalConfiguration,
                        const std::string& finalConfigurationFilename,
                        const t_inputrec*  inputrec,
                        const t_mdatoms*   mdatoms,
                        const gmx_mtop_t&  globalTop);

    //! Destructor (allows forward declaration of internal type)
    ~StatePropagatorData();

    // Allow access to state
    //! Get write access to position vector
    ArrayRefWithPadding<RVec> positionsView();
    //! Get read access to position vector
    ArrayRefWithPadding<const RVec> constPositionsView() const;
    //! Get write access to previous position vector
    ArrayRefWithPadding<RVec> previousPositionsView();
    //! Get read access to previous position vector
    ArrayRefWithPadding<const RVec> constPreviousPositionsView() const;
    //! Get write access to velocity vector
    ArrayRefWithPadding<RVec> velocitiesView();
    //! Get read access to velocity vector
    ArrayRefWithPadding<const RVec> constVelocitiesView() const;
    //! Get write access to force vector
    ForceBuffersView& forcesView();
    //! Get read access to force vector
    const ForceBuffersView& constForcesView() const;
    //! Get pointer to box
    rvec* box();
    //! Get const pointer to box
    const rvec* constBox() const;
    //! Get pointer to previous box
    rvec* previousBox();
    //! Get const pointer to previous box
    const rvec* constPreviousBox() const;
    //! Get the local number of atoms
    int localNumAtoms() const;
    //! Get the total number of atoms
    int totalNumAtoms() const;

    //! The element taking part in the simulator loop
    class Element;
    //! Get pointer to element (whose lifetime is managed by this)
    Element* element();
    //! Initial set up for the associated element
    void setup();

    //! Update the reference temperature
    void updateReferenceTemperature(ArrayRef<const real>                temperatures,
                                    ReferenceTemperatureChangeAlgorithm algorithm);
    //! Helper class handling reference temperature change
    class ReferenceTemperatureHelper;

    //! Read everything that can be stored in t_trxframe from a checkpoint file
    static void readCheckpointToTrxFrame(t_trxframe* trxFrame, ReadCheckpointData readCheckpointData);
    //! CheckpointHelper identifier
    static const std::string& checkpointID();

    //! \cond
    // (doxygen doesn't like these)
    // Classes which need access to legacy state
    friend class DomDecHelper;
    //! \endcond

private:
    //! Default constructor - only used internally
    StatePropagatorData() = default;
    //! The total number of atoms in the system
    int totalNumAtoms_;
    //! The local number of atoms
    int localNAtoms_;
    //! The position vector
    PaddedHostVector<RVec> x_;
    //! The position vector of the previous step
    PaddedHostVector<RVec> previousX_;
    //! The velocity vector
    PaddedHostVector<RVec> v_;
    //! The force vector
    ForceBuffers f_;
    //! The box matrix
    matrix box_;
    //! The box matrix of the previous step
    matrix previousBox_;
    //! The DD partitioning count
    int ddpCount_;
    //! The DD partitioning count for index_gl
    int ddpCountCgGl_;
    //! The global cg number of the local cgs
    std::vector<int> cgGl_;

    //! The global position vector
    PaddedHostVector<RVec> xGlobal_;
    //! The global position vector of the previous step
    PaddedHostVector<RVec> previousXGlobal_;
    //! The global velocity vector
    PaddedHostVector<RVec> vGlobal_;
    //! The global force vector
    PaddedHostVector<RVec> fGlobal_;

    //! The element
    std::unique_ptr<Element> element_;
    //! Instance of reference temperature helper
    std::unique_ptr<ReferenceTemperatureHelper> referenceTemperatureHelper_;

    //! Move x_ to previousX_
    void copyPosition();
    //! OMP helper to move x_ to previousX_
    void copyPosition(int start, int end);

    //! Helper function to read from / write to CheckpointData
    template<CheckpointDataOperation operation>
    void doCheckpointData(CheckpointData<operation>* checkpointData);

    // Access to legacy state
    //! Give ownership of local state resources in legacy format
    t_state* localState();
    //! Take ownership of local state resources in legacy format
    void setLocalState(t_state* state);
    /*! \brief Deep copy the local state into the provided copy and
     * return it
     *
     * In order to minimize reallocations, this function takes as a sink
     * a local state object owned by the caller, copies the current local
     * state into it, and returns the same object via a move.
     */
    std::unique_ptr<t_state> copyLocalState(std::unique_ptr<t_state> copy);
    //! Get a pointer to the global state
    t_state* globalState();
    //! Get a force pointer
    ForceBuffers* forcePointer();

    //! Whether we're doing VV and need to reset velocities after the first half step
    bool vvResetVelocities_;
    //! Velocities backup for VV
    PaddedHostVector<RVec> velocityBackup_;
    //! Function resetting the velocities
    void resetVelocities();

    //! Whether planned total number of steps was reached (used for final output only)
    bool isRegularSimulationEnd_;
    //! The signalled last step (used for final output only)
    Step lastStep_;

    // Access to ISimulator data
    //! Full simulation state (only non-nullptr on main rank).
    t_state* globalState_;
    //! Local simulation state
    t_state* localState_;
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Element for StatePropagatorData
 *
 * The `StatePropagatorData::Element` takes part in the simulator run, as it might
 * have to save a valid state at the right moment during the
 * integration. Placing the StatePropagatorData::Element correctly is the
 * duty of the simulator builder - this might be automatized later
 * if we have enough meta-data of the variables (i.e., if
 * `StatePropagatorData` knows at which time the variables currently are,
 * and can decide when a valid state (full-time step of all
 * variables) is reached. The `StatePropagatorData::Element` is also a client of
 * both the trajectory signaller and writer - it will save a
 * state for later writeout during the simulator step if it
 * knows that trajectory writing will occur later in the step,
 * and it knows how to write to file given a file pointer by
 * the `TrajectoryElement`. It is also responsible to store
 * the state for checkpointing.
 *
 */
class StatePropagatorData::Element final :
    public ISimulatorElement,
    public ITrajectoryWriterClient,
    public ITrajectorySignallerClient,
    public ICheckpointHelperClient,
    public ILastStepSignallerClient
{
public:
    //! Constructor
    Element(StatePropagatorData* statePropagatorData,
            FILE*                fplog,
            const t_commrec*     cr,
            int                  nstxout,
            int                  nstvout,
            int                  nstfout,
            int                  nstxout_compressed,
            bool                 canMoleculesBeDistributedOverPBC,
            bool                 writeFinalConfiguration,
            std::string          finalConfigurationFilename,
            const t_inputrec*    inputrec,
            const gmx_mtop_t&    globalTop);

    /*! \brief Register run function for step / time
     *
     * This needs to be called during the integration part of the simulator,
     * at the moment at which the state is at a full time step. Positioning
     * this element is the responsibility of the programmer writing the
     * integration algorithm! If the current step is a trajectory writing
     * step, StatePropagatorData will save a backup for later writeout.
     *
     * This is also the place at which the current state becomes the previous
     * state.
     *
     * \param step                 The step number
     * \param time                 The time
     * \param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    /*! \brief Backup starting velocities
     *
     * This is only needed for vv, where the first (velocity) half step is only
     * used to compute the constraint virial, but the velocities need to be reset
     * after.
     * TODO: There must be a more elegant solution to this!
     */
    void elementSetup() override;

    //! No element teardown needed
    void elementTeardown() override {}

    //! Set free energy data
    void setFreeEnergyPerturbationData(FreeEnergyPerturbationData* freeEnergyPerturbationData);

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
    //! Pointer to the associated StatePropagatorData
    StatePropagatorData* statePropagatorData_;

    //! The position writeout frequency
    const int nstxout_;
    //! The velocity writeout frequency
    const int nstvout_;
    //! The force writeout frequency
    const int nstfout_;
    //! The compressed position writeout frequency
    const int nstxout_compressed_;

    //! Pointer to keep a backup of the state for later writeout
    std::unique_ptr<t_state> localStateBackup_;
    /*! \brief Whether the contents of localStateBackup_ are logically valid
     *
     * This ensures that we don't make a second backup without consuming the
     * first. */
    bool localStateBackupValid_ = false;
    //! Step at which next writeout occurs
    Step writeOutStep_;
    //! Backup current state
    void saveState();

    //! ITrajectorySignallerClient implementation
    std::optional<SignallerCallback> registerTrajectorySignallerCallback(TrajectoryEvent event) override;

    //! ITrajectoryWriterClient implementation
    std::optional<ITrajectoryWriterCallback> registerTrajectoryWriterCallback(TrajectoryEvent event) override;

    //! ILastStepSignallerClient implementation (used for final output only)
    std::optional<SignallerCallback> registerLastStepCallback() override;

    //! Callback writing the state to file
    void write(gmx_mdoutf* outf, Step step, Time time);

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the free energy perturbation data (for trajectory writing only)
    FreeEnergyPerturbationData* freeEnergyPerturbationData_;

    //! No trajectory writer setup needed
    void trajectoryWriterSetup(gmx_mdoutf gmx_unused* outf) override {}
    //! Trajectory writer teardown - write final coordinates
    void trajectoryWriterTeardown(gmx_mdoutf* outf) override;
    //! A dummy CheckpointData - remove when we stop using the legacy trajectory writing function
    WriteCheckpointDataHolder dummyCheckpointDataHolder_;

    //! Whether planned total number of steps was reached (used for final output only)
    bool isRegularSimulationEnd_;
    //! The signalled last step (used for final output only)
    Step lastStep_;
    //! Whether system can have molecules distributed over PBC boundaries (used for final output only)
    const bool canMoleculesBeDistributedOverPBC_;
    //! Whether system has molecules self-interacting through PBC (used for final output only)
    const bool systemHasPeriodicMolecules_;
    //! The PBC type (used for final output only)
    const PbcType pbcType_;
    //! The (planned) last step - determines whether final configuration is written (used for final output only)
    const Step lastPlannedStep_;
    //! Whether final configuration was chosen in mdrun options (used for final output only)
    const bool writeFinalConfiguration_;
    //! The filename of the final configuration file (used for final output only)
    const std::string finalConfigurationFilename_;

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Handles communication.
    const t_commrec* cr_;
    //! Full system topology.
    const gmx_mtop_t& top_global_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_STATEPROPAGATORDATA_H
