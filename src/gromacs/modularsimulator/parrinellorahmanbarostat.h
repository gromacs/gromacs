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
/*! \libinternal \file
 * \brief Declares the Parrinello-Rahman barostat for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#ifndef GMX_MODULARSIMULATOR_PARRINELLORAHMANBAROSTAT_H
#define GMX_MODULARSIMULATOR_PARRINELLORAHMANBAROSTAT_H

#include "gromacs/math/vectypes.h"

#include "modularsimulatorinterfaces.h"
#include "propagator.h"

struct t_inputrec;
struct t_commrec;

namespace gmx
{
class EnergyElement;
class MDAtoms;
class StatePropagatorData;

/*! \libinternal
 * \ingroup module_modularsimulator
 * \brief Element implementing the Parrinello-Rahman barostat
 *
 * This element
 *   * integrates the Parrinello-Rahman box velocity equations,
 *   * takes a callback to the propagator to update the velocity
 *     scaling factor, and
 *   * scales the box and the positions of the system.
 */
class ParrinelloRahmanBarostat final : public ISimulatorElement, public ICheckpointHelperClient
{
public:
    //! Constructor
    ParrinelloRahmanBarostat(int                   nstpcouple,
                             int                   offset,
                             real                  couplingTimeStep,
                             Step                  initStep,
                             ArrayRef<rvec>        scalingTensor,
                             PropagatorCallbackPtr propagatorCallback,
                             StatePropagatorData*  statePropagatorData,
                             EnergyElement*        energyElement,
                             FILE*                 fplog,
                             const t_inputrec*     inputrec,
                             const MDAtoms*        mdAtoms,
                             const t_state*        globalState,
                             t_commrec*            cr,
                             bool                  isRestart);

    /*! \brief Register run function for step / time
     *
     * @param step                 The step number
     * @param time                 The time
     * @param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunctionPtr& registerRunFunction) override;

    //! Fix relative box shape
    void elementSetup() override;
    //! No element teardown needed
    void elementTeardown() override {}

    //! Getter for the box velocities
    const rvec* boxVelocities() const;

private:
    //! The frequency at which the barostat is applied
    const int nstpcouple_;
    //! If != 0, offset the step at which the barostat is applied
    const int offset_;
    //! The coupling time step - simulation time step x nstcouple_
    const real couplingTimeStep_;
    //! The first step of the simulation
    const Step initStep_;

    //! View on the velocity scaling tensor (owned by the propagator)
    ArrayRef<rvec> scalingTensor_;
    //! Callback to let propagator know that we updated lambda
    PropagatorCallbackPtr propagatorCallback_;

    //! Relative change in box before - after barostatting
    matrix mu_;
    //! Relative box shape
    tensor boxRel_;
    //! Box velocity
    tensor boxVelocity_;

    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;
    //! Pointer to the energy element
    EnergyElement* energyElement_;

    //! Integrate the PR box vector equations of motion - does not alter state
    void integrateBoxVelocityEquations(Step step);
    //! Scale box and positions
    void scaleBoxAndPositions();

    //! ICheckpointHelperClient implementation
    void writeCheckpoint(t_state* localState, t_state* globalState) override;

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Atom parameters for this domain.
    const MDAtoms* mdAtoms_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_PARRINELLORAHMANBAROSTAT_H
