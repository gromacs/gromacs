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
 * \brief Declares the propagator element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_PROPAGATOR_H
#define GMX_MODULARSIMULATOR_PROPAGATOR_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "modularsimulatorinterfaces.h"

struct gmx_wallcycle;

namespace gmx
{
class EnergyData;
class FreeEnergyPerturbationData;
class GlobalCommunicationHelper;
class LegacySimulatorData;
class MDAtoms;
class ModularSimulatorAlgorithmBuilderHelper;
class ObservablesReducer;
class StatePropagatorData;

//! \addtogroup module_modularsimulator
//! \{

//! Which velocities the thermostat scales
enum class ScaleVelocities
{
    PreStepOnly,
    PreStepAndPostStep
};

//! The different integration types we know about
enum class IntegrationStage
{
    PositionsOnly,                        //!< Moves the position vector by the given time step
    VelocitiesOnly,                       //!< Moves the velocity vector by the given time step
    LeapFrog,                             //!< Manual fusion of the previous two propagators
    VelocityVerletPositionsAndVelocities, //!< Manual position (full dt) and velocity (half dt) fusion
    ScaleVelocities,                      //!< Only scale velocities, don't propagate
    ScalePositions,                       //!< Only scale positions, don't propagate
    Count                                 //!< The number of enum entries
};

//! Sets the number of different position scaling values
enum class NumPositionScalingValues
{
    None,     //!< No position scaling (either this step or ever)
    Single,   //!< Single scaling value (either one group or all values =1)
    Multiple, //!< Multiple scaling values, need to use T-group indices
    Count     //!< The number of enum entries
};

//! Sets the number of different velocity scaling values
enum class NumVelocityScalingValues
{
    None,     //!< No velocity scaling (either this step or ever)
    Single,   //!< Single T-scaling value (either one group or all values =1)
    Multiple, //!< Multiple T-scaling values, need to use T-group indices
    Count
};

//! Describes the properties of the Parrinello-Rahman pressure scaling matrix
enum class ParrinelloRahmanVelocityScaling
{
    No,          //!< Do not apply velocity scaling (not a PR-coupling run or step)
    Diagonal,    //!< Apply velocity scaling using a diagonal matrix
    Anisotropic, //!< Apply velocity scaling using a matrix with off-diagonal elements
    Count
};

/*! \internal
 * \brief Propagator element
 *
 * The propagator element can, through templating, cover the different
 * propagation types used in NVE MD. The combination of templating, static
 * functions, and having only the inner-most operations in the static
 * functions allows to have performance comparable to fused update elements
 * while keeping easily re-orderable single instructions.
 *
 * \tparam integrationStep  The integration types
 */
template<IntegrationStage integrationStage>
class Propagator final : public ISimulatorElement
{
public:
    //! Constructor
    Propagator(double               timestep,
               StatePropagatorData* statePropagatorData,
               const MDAtoms*       mdAtoms,
               gmx_wallcycle*       wcycle);

    /*! \brief Register run function for step / time
     *
     * This function will determine the required flavor of the run function to be registered
     * for the current step. In case of the pure scaling integrator stage, it might also skip
     * the function registration if no scaling is needed.
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

    //! Set the number of velocity scaling variables
    void setNumVelocityScalingVariables(int numVelocityScalingVariables, ScaleVelocities scaleVelocities);
    //! Set the number of position scaling variables
    void setNumPositionScalingVariables(int numPositionScalingVariables);
    //! Get view on the scaling vector applied to start of step velocities
    ArrayRef<real> viewOnStartVelocityScaling();
    //! Get view on the scaling vector applied to end of step velocities
    ArrayRef<real> viewOnEndVelocityScaling();
    //! Get view on the scaling vector applied to the positions
    ArrayRef<real> viewOnPositionScaling();
    //! Get velocity scaling callback
    PropagatorCallback velocityScalingCallback();
    //! Get position scaling callback
    PropagatorCallback positionScalingCallback();

    //! Get view on the full PR scaling matrix
    ArrayRef<rvec> viewOnPRScalingMatrix();
    //! Get PR scaling callback
    PropagatorCallback prScalingCallback();

    /*! \brief Factory method implementation
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper   Pointer to the \c GlobalCommunicationHelper object
     * \param observablesReducer          Pointer to the \c ObservablesReducer object
     * \param propagatorTag  The name of the propagator to simplify connection
     * \param timestep  The time step the propagator uses
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement* getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                    StatePropagatorData*        statePropagatorData,
                                                    EnergyData*                 energyData,
                                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                    GlobalCommunicationHelper* globalCommunicationHelper,
                                                    ObservablesReducer*        observablesReducer,
                                                    const PropagatorTag&       propagatorTag,
                                                    TimeStep                   timestep);

    /*! \brief Factory method implementation
     *
     * Version without time step for pure scaling elements
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper   Pointer to the \c GlobalCommunicationHelper object
     * \param observablesReducer          Pointer to the \c ObservablesReducer object
     * \param propagatorTag  The name of the propagator to simplify connection
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement* getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                    StatePropagatorData*        statePropagatorData,
                                                    EnergyData*                 energyData,
                                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                    GlobalCommunicationHelper* globalCommunicationHelper,
                                                    ObservablesReducer*        observablesReducer,
                                                    const PropagatorTag&       propagatorTag);

private:
    //! The actual propagation
    template<NumVelocityScalingValues        numStartVelocityScalingValues,
             ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
             NumVelocityScalingValues        numEndVelocityScalingValues,
             NumPositionScalingValues        numPositionScalingValues>
    void run();

    //! The time step
    const real timestep_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;

    //! Whether we're doing single-value velocity scaling (velocities at start of step)
    bool doSingleStartVelocityScaling_;
    //! Whether we're doing group-wise velocity scaling (velocities at start of step)
    bool doGroupStartVelocityScaling_;
    //! Whether we're doing single-value velocity scaling (velocities at end of step)
    bool doSingleEndVelocityScaling_;
    //! Wether we're doing group-wise velocity scaling (velocities at end of step)
    bool doGroupEndVelocityScaling_;
    //! Whether we're doing single-value position scaling
    bool doSinglePositionScaling_;
    //! Whether we're doing group-wise position scaling
    bool doGroupPositionScaling_;
    //! The vector of velocity scaling values
    std::vector<real> startVelocityScaling_;
    //! The vector of velocity scaling values
    std::vector<real> endVelocityScaling_;
    //! The vector of position scaling values
    std::vector<real> positionScaling_;
    //! The next velocity scaling step
    Step scalingStepVelocity_;
    //! The next position scaling step
    Step scalingStepPosition_;

    //! The diagonal of the PR scaling matrix
    rvec diagPR_;
    //! The full PR scaling matrix
    matrix matrixPR_;
    //! The next PR scaling step
    Step scalingStepPR_;

    // Access to ISimulator data
    //! Atom parameters for this domain.
    const MDAtoms* mdAtoms_;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wcycle_;
};

//! \}
} // namespace gmx

#endif // GMX_MODULARSIMULATOR_PROPAGATOR_H
