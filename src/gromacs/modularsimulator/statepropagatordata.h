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
 * \brief Declares the state for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#ifndef GMX_MODULARSIMULATOR_STATEPROPAGATORDATA_H
#define GMX_MODULARSIMULATOR_STATEPROPAGATORDATA_H

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"

#include "modularsimulatorinterfaces.h"
#include "topologyholder.h"

struct gmx_mdoutf;
struct t_commrec;
struct t_inputrec;
class t_state;
struct t_mdatoms;

namespace gmx
{
enum class ConstraintVariable;
class FreeEnergyPerturbationElement;

/*! \libinternal
 * \ingroup module_modularsimulator
 * \brief StatePropagatorData and associated data
 *
 * The `StatePropagatorData` contains a little more than the pure
 * statistical-physical micro state, namely the positions,
 * velocities, forces, and box matrix, as well as a backup of
 * the positions and box of the last time step. While it takes
 * part in the simulator loop to be able to backup positions /
 * boxes and save the current state if needed, it's main purpose
 * is to offer access to its data via getter methods. All elements
 * reading or writing to this data need a pointer to the
 * `StatePropagatorData` and need to request their data explicitly. This
 * will later simplify the understanding of data dependencies
 * between elements.
 *
 * The `StatePropagatorData` takes part in the simulator run, as it might
 * have to save a valid state at the right moment during the
 * integration. Placing the StatePropagatorData correctly is for now the
 * duty of the simulator builder - this might be automatized later
 * if we have enough meta-data of the variables (i.e., if
 * `StatePropagatorData` knows at which time the variables currently are,
 * and can decide when a valid state (full-time step of all
 * variables) is reached. The `StatePropagatorData` is also a client of
 * both the trajectory signaller and writer - it will save a
 * state for later writeout during the simulator step if it
 * knows that trajectory writing will occur later in the step,
 * and it knows how to write to file given a file pointer by
 * the `TrajectoryElement`.
 *
 * Note that the `StatePropagatorData` can be converted to and from the
 * legacy `t_state` object. This is useful when dealing with
 * functionality which has not yet been adapted to use the new
 * data approach - of the elements currently implemented, only
 * domain decomposition, PME load balancing, and the initial
 * constraining are using this.
 */
class StatePropagatorData final :
    public ISimulatorElement,
    public ITrajectoryWriterClient,
    public ITrajectorySignallerClient,
    public ICheckpointHelperClient,
    public ILastStepSignallerClient
{
public:
    //! Constructor
    StatePropagatorData(int                            numAtoms,
                        FILE*                          fplog,
                        const t_commrec*               cr,
                        t_state*                       globalState,
                        int                            nstxout,
                        int                            nstvout,
                        int                            nstfout,
                        int                            nstxout_compressed,
                        bool                           useGPU,
                        FreeEnergyPerturbationElement* freeEnergyPerturbationElement,
                        const TopologyHolder*          topologyHolder,
                        bool                           canMoleculesBeDistributedOverPBC,
                        bool                           writeFinalConfiguration,
                        std::string                    finalConfigurationFilename,
                        const t_inputrec*              inputrec,
                        const t_mdatoms*               mdatoms);

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
    ArrayRefWithPadding<RVec> forcesView();
    //! Get read access to force vector
    ArrayRefWithPadding<const RVec> constForcesView() const;
    //! Get pointer to box
    rvec* box();
    //! Get const pointer to box
    const rvec* constBox();
    //! Get pointer to previous box
    rvec* previousBox();
    //! Get const pointer to previous box
    const rvec* constPreviousBox();
    //! Get the local number of atoms
    int localNumAtoms();

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
     * @param step                 The step number
     * @param time                 The time
     * @param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunctionPtr& registerRunFunction) override;

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

    //! @cond
    // (doxygen doesn't like these)
    // Classes which need access to legacy state
    friend class DomDecHelper;
    //! @endcond

private:
    //! The total number of atoms in the system
    int totalNumAtoms_;
    //! The position writeout frequency
    int nstxout_;
    //! The velocity writeout frequency
    int nstvout_;
    //! The force writeout frequency
    int nstfout_;
    //! The compressed position writeout frequency
    int nstxout_compressed_;

    //! The local number of atoms
    int localNAtoms_;
    //! The position vector
    PaddedHostVector<RVec> x_;
    //! The position vector of the previous step
    PaddedHostVector<RVec> previousX_;
    //! The velocity vector
    PaddedHostVector<RVec> v_;
    //! The force vector
    PaddedHostVector<RVec> f_;
    //! The box matrix
    matrix box_;
    //! The box matrix of the previous step
    matrix previousBox_;
    //! The DD partitioning count for legacy t_state compatibility
    int ddpCount_;

    //! Move x_ to previousX_
    void copyPosition();
    //! OMP helper to move x_ to previousX_
    void copyPosition(int start, int end);

    // Access to legacy state
    //! Get a deep copy of the current state in legacy format
    std::unique_ptr<t_state> localState();
    //! Update the current state with a state in legacy format
    void setLocalState(std::unique_ptr<t_state> state);
    //! Get a pointer to the global state
    t_state* globalState();
    //! Get a force pointer
    PaddedHostVector<gmx::RVec>* forcePointer();

    //! Pointer to keep a backup of the state for later writeout
    std::unique_ptr<t_state> localStateBackup_;
    //! Step at which next writeout occurs
    Step writeOutStep_;
    //! Backup current state
    void saveState();

    //! ITrajectorySignallerClient implementation
    SignallerCallbackPtr registerTrajectorySignallerCallback(TrajectoryEvent event) override;

    //! ITrajectoryWriterClient implementation
    ITrajectoryWriterCallbackPtr registerTrajectoryWriterCallback(TrajectoryEvent event) override;

    //! ICheckpointHelperClient implementation
    void writeCheckpoint(t_state* localState, t_state* globalState) override;

    //! ILastStepSignallerClient implementation (used for final output only)
    SignallerCallbackPtr registerLastStepCallback() override;

    //! Callback writing the state to file
    void write(gmx_mdoutf* outf, Step step, Time time);

    //! Whether we're doing VV and need to reset velocities after the first half step
    bool vvResetVelocities_;
    //! Velocities backup for VV
    PaddedHostVector<RVec> velocityBackup_;
    //! Function resetting the velocities
    void resetVelocities();

    //! Pointer to the free energy perturbation element (for trajectory writing only)
    FreeEnergyPerturbationElement* freeEnergyPerturbationElement_;

    //! Whether planned total number of steps was reached (used for final output only)
    bool isRegularSimulationEnd_;
    //! The signalled last step (used for final output only)
    Step lastStep_;

    //! Whether system can have molecules distributed over PBC boundaries (used for final output only)
    const bool canMoleculesBeDistributedOverPBC_;
    //! Whether system has molecules self-interacting through PBC (used for final output only)
    const bool systemHasPeriodicMolecules_;
    //! The PBC type (used for final output only)
    const int pbcType_;
    //! Pointer to the topology (used for final output only)
    const TopologyHolder* topologyHolder_;
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
    //! Full simulation state (only non-nullptr on master rank).
    t_state* globalState_;

    //! No trajectory writer setup needed
    void trajectoryWriterSetup(gmx_mdoutf gmx_unused* outf) override {}
    //! Trajectory writer teardown - write final coordinates
    void trajectoryWriterTeardown(gmx_mdoutf* outf) override;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_STATEPROPAGATORDATA_H
