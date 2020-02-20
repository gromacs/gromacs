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
 * \brief Declares the shell / flex constraints element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#ifndef GMX_MODULARSIMULATOR_SHELLFCELEMENT_H
#define GMX_MODULARSIMULATOR_SHELLFCELEMENT_H

#include <array>

#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/real.h"

#include "modularsimulatorinterfaces.h"
#include "topologyholder.h"

struct gmx_enfrot;
struct gmx_shellfc_t;
struct gmx_wallcycle;
struct pull_t;
struct t_fcdata;
struct t_nrnb;

namespace gmx
{
class Awh;
class EnergyElement;
class FreeEnergyPerturbationElement;
class ImdSession;
class MDAtoms;
class MdrunScheduleWorkload;
class StatePropagatorData;

/*! \libinternal
 * \ingroup module_modularsimulator
 * \brief Shell & flex constraints element
 *
 * The ShellFCElement manages the call to relax_shell_flexcon(...)
 */
class ShellFCElement final :
    public ISimulatorElement,
    public ITopologyHolderClient,
    public INeighborSearchSignallerClient,
    public IEnergySignallerClient
{
public:
    //! Constructor
    ShellFCElement(StatePropagatorData*           statePropagatorData,
                   EnergyElement*                 energyElement,
                   FreeEnergyPerturbationElement* freeEnergyPerturbationElement,
                   bool                           isVerbose,
                   bool                           isDynamicBox,
                   FILE*                          fplog,
                   const t_commrec*               cr,
                   const t_inputrec*              inputrec,
                   const MDAtoms*                 mdAtoms,
                   t_nrnb*                        nrnb,
                   t_forcerec*                    fr,
                   t_fcdata*                      fcd,
                   gmx_wallcycle*                 wcycle,
                   MdrunScheduleWorkload*         runScheduleWork,
                   gmx_vsite_t*                   vsite,
                   ImdSession*                    imdSession,
                   pull_t*                        pull_work,
                   Constraints*                   constr,
                   const gmx_mtop_t*              globalTopology,
                   gmx_enfrot*                    enforcedRotation);

    /*! \brief Register shell / flex constraint calculation for step / time
     *
     * @param step                 The step number
     * @param time                 The time
     * @param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunctionPtr& registerRunFunction) override;

    //! Check that we got the local topology
    void elementSetup() override;
    //! Print some final output
    void elementTeardown() override;

    //! Whether either shells or flexible constraints are used
    static bool doShellsOrFlexConstraints(const gmx_mtop_t& mtop, int nflexcon);

private:
    //! ITopologyHolderClient implementation
    void setTopology(const gmx_localtop_t* top) override;
    //! INeighborSearchSignallerClient implementation
    SignallerCallbackPtr registerNSCallback() override;
    //! IEnergySignallerClient implementation
    SignallerCallbackPtr registerEnergyCallback(EnergySignallerEvent event) override;
    //! The actual do_force call
    void run(Step step, Time time, bool isNSStep, unsigned int flags);

    //! The shell / FC helper struct
    gmx_shellfc_t* shellfc_;

    //! The next NS step
    Step nextNSStep_;
    //! The next energy calculation step
    Step nextEnergyCalculationStep_;
    //! The next energy calculation step
    Step nextVirialCalculationStep_;
    //! The next free energy calculation step
    Step nextFreeEnergyCalculationStep_;

    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;
    //! Pointer to the energy element
    EnergyElement* energyElement_;
    //! The free energy perturbation element
    FreeEnergyPerturbationElement* freeEnergyPerturbationElement_;

    //! The local topology - updated by Topology via Client system
    const gmx_localtop_t* localTopology_;

    //! Whether we're having a dynamic box
    const bool isDynamicBox_;
    //! Whether we're being verbose
    const bool isVerbose_;
    //! The number of shell relaxation steps we did
    Step nSteps_;

    //! DD / DLB helper object
    const DDBalanceRegionHandler ddBalanceRegionHandler_;

    /* \brief The FEP lambda vector
     *
     * Used if FEP is off, since do_force
     * requires lambda to be allocated anyway
     */
    std::array<real, efptNR> lambda_;

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Handles communication.
    const t_commrec* cr_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Atom parameters for this domain.
    const MDAtoms* mdAtoms_;
    //! Manages flop accounting.
    t_nrnb* nrnb_;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wcycle_;
    //! Parameters for force calculations.
    t_forcerec* fr_;
    //! Handles virtual sites.
    gmx_vsite_t* vsite_;
    //! The Interactive Molecular Dynamics session.
    ImdSession* imdSession_;
    //! The pull work object.
    pull_t* pull_work_;
    //! Helper struct for force calculations.
    t_fcdata* fcd_;
    //! Schedule of work for each MD step for this task.
    MdrunScheduleWorkload* runScheduleWork_;
    //! Handles constraints.
    Constraints* constr_;
    //! Handles enforced rotation.
    gmx_enfrot* enforcedRotation_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_SHELLFCELEMENT_H
