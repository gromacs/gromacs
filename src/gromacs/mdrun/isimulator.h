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
struct gmx_ekindata_t;
struct gmx_enerdata_t;
struct gmx_enfrot;
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
struct MdModulesNotifier;
class ImdSession;
class MDLogger;
class MDAtoms;
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
                        const MdModulesNotifier&            mdModulesNotifier,
                        t_inputrec*                         inputrec,
                        ImdSession*                         imdSession,
                        pull_t*                             pull_work,
                        t_swap*                             swap,
                        gmx_mtop_t*                         top_global,
                        t_state*                            state_global,
                        ObservablesHistory*                 observablesHistory,
                        MDAtoms*                            mdAtoms,
                        t_nrnb*                             nrnb,
                        gmx_wallcycle*                      wcycle,
                        t_forcerec*                         fr,
                        gmx_enerdata_t*                     enerd,
                        gmx_ekindata_t*                     ekind,
                        MdrunScheduleWorkload*              runScheduleWork,
                        const ReplicaExchangeParameters&    replExParams,
                        gmx_membed_t*                       membed,
                        gmx_walltime_accounting*            walltime_accounting,
                        std::unique_ptr<StopHandlerBuilder> stopHandlerBuilder,
                        bool                                doRerun) :
        fplog(fplog),
        cr(cr),
        ms(ms),
        mdlog(mdlog),
        nfile(nfile),
        fnm(fnm),
        oenv(oenv),
        mdrunOptions(mdrunOptions),
        startingBehavior(startingBehavior),
        vsite(vsite),
        constr(constr),
        enforcedRotation(enforcedRotation),
        deform(deform),
        outputProvider(outputProvider),
        mdModulesNotifier(mdModulesNotifier),
        inputrec(inputrec),
        imdSession(imdSession),
        pull_work(pull_work),
        swap(swap),
        top_global(top_global),
        state_global(state_global),
        observablesHistory(observablesHistory),
        mdAtoms(mdAtoms),
        nrnb(nrnb),
        wcycle(wcycle),
        fr(fr),
        enerd(enerd),
        ekind(ekind),
        runScheduleWork(runScheduleWork),
        replExParams(replExParams),
        membed(membed),
        walltime_accounting(walltime_accounting),
        stopHandlerBuilder(std::move(stopHandlerBuilder)),
        doRerun(doRerun)
    {
    }

    //! Handles logging.
    FILE* fplog;
    //! Handles communication.
    t_commrec* cr;
    //! Coordinates multi-simulations.
    const gmx_multisim_t* ms;
    //! Handles logging.
    const MDLogger& mdlog;
    //! Count of input file options.
    int nfile;
    //! Content of input file options.
    const t_filenm* fnm;
    //! Handles writing text output.
    const gmx_output_env_t* oenv;
    //! Contains command-line options to mdrun.
    const MdrunOptions& mdrunOptions;
    //! Whether the simulation will start afresh, or restart with/without appending.
    const StartingBehavior startingBehavior;
    //! Handles virtual sites.
    VirtualSitesHandler* vsite;
    //! Handles constraints.
    Constraints* constr;
    //! Handles enforced rotation.
    gmx_enfrot* enforcedRotation;
    //! Handles box deformation.
    BoxDeformation* deform;
    //! Handles writing output files.
    IMDOutputProvider* outputProvider;
    //! Handles notifications to MdModules for checkpoint writing
    const MdModulesNotifier& mdModulesNotifier;
    //! Contains user input mdp options.
    t_inputrec* inputrec;
    //! The Interactive Molecular Dynamics session.
    ImdSession* imdSession;
    //! The pull work object.
    pull_t* pull_work;
    //! The coordinate-swapping session.
    t_swap* swap;
    //! Full system topology.
    const gmx_mtop_t* top_global;
    //! Full simulation state (only non-nullptr on master rank).
    t_state* state_global;
    //! History of simulation observables.
    ObservablesHistory* observablesHistory;
    //! Atom parameters for this domain.
    MDAtoms* mdAtoms;
    //! Manages flop accounting.
    t_nrnb* nrnb;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wcycle;
    //! Parameters for force calculations.
    t_forcerec* fr;
    //! Data for energy output.
    gmx_enerdata_t* enerd;
    //! Kinetic energy data.
    gmx_ekindata_t* ekind;
    //! Schedule of work for each MD step for this task.
    MdrunScheduleWorkload* runScheduleWork;
    //! Parameters for replica exchange algorihtms.
    const ReplicaExchangeParameters& replExParams;
    //! Parameters for membrane embedding.
    gmx_membed_t* membed;
    //! Manages wall time accounting.
    gmx_walltime_accounting* walltime_accounting;
    //! Registers stop conditions
    std::unique_ptr<StopHandlerBuilder> stopHandlerBuilder;
    //! Whether we're doing a rerun.
    bool doRerun;
};

} // namespace gmx

#endif // GMX_MDRUN_ISIMULATOR_H
