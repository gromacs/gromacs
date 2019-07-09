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
/*! \libinternal
 * \brief Declares the propagator element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#ifndef GMX_MODULARSIMULATOR_PROPAGATOR_H
#define GMX_MODULARSIMULATOR_PROPAGATOR_H

#include "gromacs/utility/real.h"

#include "modularsimulatorinterfaces.h"

struct gmx_wallcycle;

namespace gmx
{
class MDAtoms;
class StatePropagatorData;

//! \addtogroup module_modularsimulator
//! \{

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

/*! \libinternal
 * \brief Propagator element
 *
 * The propagator element can, through templating, cover the different
 * propagation types used in NVE MD. The combination of templating, static
 * functions, and having only the inner-most operations in the static
 * functions allows to have performance comparable to fused update elements
 * while keeping easily re-orderable single instructions.
 *
 * \todo: Get rid of updateVelocities2() once we don't require identical
 *        reproduction of do_md() results.
 *
 * @tparam algorithm  The integration types
 */
template <IntegrationStep algorithm>
class Propagator final :
    public       ISimulatorElement
{
    public:
        //! Constructor
        Propagator(
            double               timestep,
            StatePropagatorData *statePropagatorData,
            const MDAtoms       *mdAtoms,
            gmx_wallcycle       *wcycle);

        /*! \brief Register run function for step / time
         *
         * @param step                 The step number
         * @param time                 The time
         * @param registerRunFunction  Function allowing to register a run function
         */
        void scheduleTask(
            Step step, Time time,
            const RegisterRunFunctionPtr &registerRunFunction) override;

        //! No element setup needed
        void elementSetup() override {}
        //! No element teardown needed
        void elementTeardown() override {}

    private:
        //! The actual propagation
        void run();

        //! The time step
        const real timestep_;

        //! Pointer to the micro state
        StatePropagatorData *statePropagatorData_;

        // Access to ISimulator data
        //! Atom parameters for this domain.
        const MDAtoms *mdAtoms_;
        //! Manages wall cycle accounting.
        gmx_wallcycle *wcycle_;

};

//! \}
}      // namespace gmx

#endif // GMX_MODULARSIMULATOR_PROPAGATOR_H
