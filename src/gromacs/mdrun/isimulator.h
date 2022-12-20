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
/*! \internal
 * \brief Declares the general simulator interface
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_ISIMULATOR_H
#define GMX_MDRUN_ISIMULATOR_H

#include "gromacs/mdlib/stophandler.h"

class energyhistory_t;
class gmx_ekindata_t;
struct gmx_enerdata_t;
struct gmx_enfrot;
struct gmx_localtop_t;
struct gmx_mtop_t;
struct gmx_membed_t;
struct gmx_multisim_t;
struct gmx_output_env_t;
struct gmx_wallcycle;
struct gmx_walltime_accounting;
struct ObservablesHistory;
struct pull_t;
struct ReplicaExchangeParameters;
struct t_commrec;
struct t_forcerec;
struct t_filenm;
struct t_inputrec;
struct t_nrnb;
struct t_swap;
class t_state;

namespace gmx
{
enum class StartingBehavior;
class BoxDeformation;
class Constraints;
class MdrunScheduleWorkload;
class IMDOutputProvider;
struct MDModulesNotifiers;
class ImdSession;
class MDLogger;
class MDAtoms;
class ObservablesReducerBuilder;
class StopHandlerBuilder;
struct MdrunOptions;
class VirtualSitesHandler;

/*! \internal
 * \brief The Simulator interface
 *
 * This is the general interface for any type of simulation type
 * run with GROMACS. This allows having a builder return different
 * Simulator objects based on user input.
 */
class ISimulator
{
public:
    /*! \brief The simulation run
     *
     * This will be called by the owner of the simulator object. To be redefined
     * by the child classes. This function is expected to run the simulation.
     */
    virtual void run() = 0;
    //! Standard destructor
    virtual ~ISimulator() = default;
};

/*! \internal
 * \brief The legacy simulator data
 *
 * This contains the data passed into the GROMACS simulators from
 * the Mdrunner object.
 */
class LegacySimulatorData
{
public:
    //! The constructor
    LegacySimulatorData(FILE*                               fplog,
                        t_commrec*                          cr,
                        const gmx_multisim_t*               ms,
                        const MDLogger&                     mdlog,
                        int                                 nfile,
                        const t_filenm*                     fnm,
                        const gmx_output_env_t*             oenv,
                        const MdrunOptions&                 mdrunOptions,
                        StartingBehavior                    startingBehavior,
                        VirtualSitesHandler*                vsite,
                        Constraints*                        constr,
                        gmx_enfrot*                         enforcedRotation,
                        BoxDeformation*                     deform,
                        IMDOutputProvider*                  outputProvider,
                        const MDModulesNotifiers&           mdModulesNotifiers,
                        t_inputrec*                         inputrec,
                        ImdSession*                         imdSession,
                        pull_t*                             pull_work,
                        t_swap*                             swap,
                        const gmx_mtop_t&                   top_global,
                        gmx_localtop_t*                     top,
                        t_state*                            state_global,
                        t_state*                            state,
                        ObservablesHistory*                 observablesHistory,
                        MDAtoms*                            mdAtoms,
                        t_nrnb*                             nrnb,
                        gmx_wallcycle*                      wcycle,
                        t_forcerec*                         fr,
                        gmx_enerdata_t*                     enerd,
                        ObservablesReducerBuilder*          observablesReducerBuilder,
                        gmx_ekindata_t*                     ekind,
                        MdrunScheduleWorkload*              runScheduleWork,
                        const ReplicaExchangeParameters&    replExParams,
                        gmx_membed_t*                       membed,
                        gmx_walltime_accounting*            walltime_accounting,
                        std::unique_ptr<StopHandlerBuilder> stopHandlerBuilder,
                        bool                                doRerun) :
        fpLog_(fplog),
        cr_(cr),
        ms_(ms),
        mdLog_(mdlog),
        nFile_(nfile),
        fnm_(fnm),
        oenv_(oenv),
        mdrunOptions_(mdrunOptions),
        startingBehavior_(startingBehavior),
        virtualSites_(vsite),
        constr_(constr),
        enforcedRotation_(enforcedRotation),
        deform_(deform),
        outputProvider_(outputProvider),
        mdModulesNotifiers_(mdModulesNotifiers),
        inputRec_(inputrec),
        imdSession_(imdSession),
        pullWork_(pull_work),
        swap_(swap),
        topGlobal_(top_global),
        top_(top),
        stateGlobal_(state_global),
        state_(state),
        observablesHistory_(observablesHistory),
        mdAtoms_(mdAtoms),
        nrnb_(nrnb),
        wallCycleCounters_(wcycle),
        fr_(fr),
        enerd_(enerd),
        observablesReducerBuilder_(observablesReducerBuilder),
        ekind_(ekind),
        runScheduleWork_(runScheduleWork),
        replExParams_(replExParams),
        membed_(membed),
        wallTimeAccounting_(walltime_accounting),
        stopHandlerBuilder_(std::move(stopHandlerBuilder)),
        doRerun_(doRerun)
    {
    }

    //! Handles logging.
    FILE* fpLog_;
    //! Handles communication.
    t_commrec* cr_;
    //! Coordinates multi-simulations.
    const gmx_multisim_t* ms_;
    //! Handles logging.
    const MDLogger& mdLog_;
    //! Count of input file options.
    int nFile_;
    //! Content of input file options.
    const t_filenm* fnm_;
    //! Handles writing text output.
    const gmx_output_env_t* oenv_;
    //! Contains command-line options to mdrun.
    const MdrunOptions& mdrunOptions_;
    //! Whether the simulation will start afresh, or restart with/without appending.
    const StartingBehavior startingBehavior_;
    //! Handles virtual sites.
    VirtualSitesHandler* virtualSites_;
    //! Handles constraints.
    Constraints* constr_;
    //! Handles enforced rotation.
    gmx_enfrot* enforcedRotation_;
    //! Handles box deformation.
    BoxDeformation* deform_;
    //! Handles writing output files.
    IMDOutputProvider* outputProvider_;
    //! Handles notifications to MDModules for checkpoint writing
    const MDModulesNotifiers& mdModulesNotifiers_;
    //! Contains user input mdp options. Note: The const-ness is casted away in a few instances, see #3854.
    const t_inputrec* inputRec_;
    //! The Interactive Molecular Dynamics session.
    ImdSession* imdSession_;
    //! The pull work object.
    pull_t* pullWork_;
    //! The coordinate-swapping session.
    t_swap* swap_;
    //! Full system topology.
    const gmx_mtop_t& topGlobal_;
    //! Handle to local simulation topology.
    gmx_localtop_t* top_;
    //! Full simulation state (only non-nullptr on main rank).
    t_state* stateGlobal_;
    //! Handle to local state of the simulation.
    t_state* state_;
    //! History of simulation observables.
    ObservablesHistory* observablesHistory_;
    //! Atom parameters for this domain.
    MDAtoms* mdAtoms_;
    //! Manages flop accounting.
    t_nrnb* nrnb_;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wallCycleCounters_;
    //! Parameters for force calculations.
    t_forcerec* fr_;
    //! Data for energy output.
    gmx_enerdata_t* enerd_;
    //! Builder for coordinator of reduction for observables
    ObservablesReducerBuilder* observablesReducerBuilder_;
    //! Kinetic energy data.
    gmx_ekindata_t* ekind_;
    //! Schedule of work for each MD step for this task.
    MdrunScheduleWorkload* runScheduleWork_;
    //! Parameters for replica exchange algorihtms.
    const ReplicaExchangeParameters& replExParams_;
    //! Parameters for membrane embedding.
    gmx_membed_t* membed_;
    //! Manages wall time accounting.
    gmx_walltime_accounting* wallTimeAccounting_;
    //! Registers stop conditions
    std::unique_ptr<StopHandlerBuilder> stopHandlerBuilder_;
    //! Whether we're doing a rerun.
    bool doRerun_;
};

} // namespace gmx

#endif // GMX_MDRUN_ISIMULATOR_H
