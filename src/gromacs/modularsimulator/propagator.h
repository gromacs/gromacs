/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
class StatePropagatorData;

//! \addtogroup module_modularsimulator
//! \{

//! Whether built propagator should be registered with thermostat
enum class RegisterWithThermostat
{
    True,
    False
};
//! Whether built propagator should be registered with barostat
enum class RegisterWithBarostat
{
    True,
    False
};

/*! \brief The different integration types we know about
 *
 * PositionsOnly:
 *   Moves the position vector by the given time step
 * VelocitiesOnly:
 *   Moves the velocity vector by the given time step
 * LeapFrog:
 *   Is a manual fusion of the previous two propagators
 * VelocityVerletPositionsAndVelocities:
 *   Is a manual fusion of VelocitiesOnly and PositionsOnly,
 *   where VelocitiesOnly is only propagated by half the
 *   time step of the positions.
 */
enum class IntegrationStep
{
    PositionsOnly,
    VelocitiesOnly,
    LeapFrog,
    VelocityVerletPositionsAndVelocities,
    Count
};

//! Sets the number of different velocity scaling values
enum class NumVelocityScalingValues
{
    None,     //!< No velocity scaling (either this step or ever)
    Single,   //!< Single T-scaling value (either one group or all values =1)
    Multiple, //!< Multiple T-scaling values, need to use T-group indices
    Count
};

//! Sets the type of Parrinello-Rahman pressure scaling
enum class ParrinelloRahmanVelocityScaling
{
    No,       //!< Do not apply velocity scaling (not a PR-coupling run or step)
    Diagonal, //!< Apply velocity scaling using a diagonal matrix
    Full,     //!< Apply velocity scaling using a full matrix
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
 * \tparam algorithm  The integration types
 */
template<IntegrationStep algorithm>
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
    void setNumVelocityScalingVariables(int numVelocityScalingVariables);
    //! Get view on the velocity scaling vector
    ArrayRef<real> viewOnVelocityScaling();
    //! Get velocity scaling callback
    PropagatorCallback velocityScalingCallback();

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
     * \param globalCommunicationHelper  Pointer to the \c GlobalCommunicationHelper object
     * \param timestep  The time step the propagator uses
     * \param registerWithThermostat  Whether this propagator should be registered with the thermostat
     * \param registerWithBarostat  Whether this propagator should be registered with the barostat
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement* getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                    StatePropagatorData*        statePropagatorData,
                                                    EnergyData*                 energyData,
                                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                    GlobalCommunicationHelper* globalCommunicationHelper,
                                                    double                     timestep,
                                                    RegisterWithThermostat registerWithThermostat,
                                                    RegisterWithBarostat   registerWithBarostat);

private:
    //! The actual propagation
    template<NumVelocityScalingValues numVelocityScalingValues, ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling>
    void run();

    //! The time step
    const real timestep_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;

    //! Whether we're doing single-value velocity scaling
    bool doSingleVelocityScaling_;
    //! Wether we're doing group-wise velocity scaling
    bool doGroupVelocityScaling_;
    //! The vector of velocity scaling values
    std::vector<real> velocityScaling_;
    //! The next velocity scaling step
    Step scalingStepVelocity_;

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
