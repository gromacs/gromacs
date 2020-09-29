/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019-2020, by the GROMACS development team, led by
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
 * \brief Declares the simulator builder for mdrun
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_SIMULATORBUILDER_H
#define GMX_MDRUN_SIMULATORBUILDER_H

#include <memory>


class energyhistory_t;
struct gmx_ekindata_t;
struct gmx_enerdata_t;
struct gmx_enfrot;
struct gmx_mtop_t;
struct gmx_multisim_t;
struct gmx_output_env_t;
struct gmx_wallcycle;
struct gmx_walltime_accounting;
struct ObservablesHistory;
struct pull_t;
struct ReplicaExchangeParameters;
struct t_commrec;
struct t_filenm;
struct t_forcerec;
struct t_inputrec;
struct t_nrnb;
class t_state;
struct t_swap;

namespace gmx
{
class BoxDeformation;
class Constraints;
class IMDOutputProvider;
class ImdSession;
class ISimulator;
class MdrunScheduleWorkload;
class MembedHolder;
class MDAtoms;
class MDLogger;
struct MdModulesNotifier;
struct MdrunOptions;
class ReadCheckpointDataHolder;
enum class StartingBehavior;
class StopHandlerBuilder;
class VirtualSitesHandler;

/*! \brief
 * Simulation configuation settings.
 */
struct SimulatorConfig
{
public:
    //! Build from settings for this simulation.
    SimulatorConfig(const MdrunOptions&    mdrunOptions,
                    StartingBehavior       startingBehavior,
                    MdrunScheduleWorkload* runScheduleWork) :
        mdrunOptions_(mdrunOptions),
        startingBehavior_(startingBehavior),
        runScheduleWork_(runScheduleWork)
    {
    }
    // TODO: Specify copy and move semantics.

    //! Handle to user options.
    const MdrunOptions& mdrunOptions_;
    //! How are we starting the simulation.
    StartingBehavior startingBehavior_;
    //! How are we scheduling the tasks for this simulation.
    MdrunScheduleWorkload* runScheduleWork_;
};


/*! \brief
 * Data for a specific simulation state.
 *
 * \todo Think of a better name and annoy people that forget
 *       to add documentation for their code.
 */
struct SimulatorStateData
{
    //! Build collection of current state data.
    SimulatorStateData(t_state*            globalState,
                       ObservablesHistory* observablesHistory,
                       gmx_enerdata_t*     enerdata,
                       gmx_ekindata_t*     ekindata) :
        globalState_p(globalState),
        observablesHistory_p(observablesHistory),
        enerdata_p(enerdata),
        ekindata_p(ekindata)
    {
    }

    //! Can perform copy of current state.
    SimulatorStateData(const SimulatorStateData& simulatorStateData) = default;

    //! Handle to global state of the simulation.
    t_state* globalState_p;
    //! Handle to current simulation history.
    ObservablesHistory* observablesHistory_p;
    //! Handle to collected data for energy groups.
    gmx_enerdata_t* enerdata_p;
    //! Handle to collected data for kinectic energy.
    gmx_ekindata_t* ekindata_p;
};

/*! \brief
 * Collection of environmental information for a simulation.
 *
 * \todo Fix doxygen checks.
 */
class SimulatorEnv
{
public:
    //! Build from current simulation environment.
    SimulatorEnv(FILE*             fplog,
                 t_commrec*        commRec,
                 gmx_multisim_t*   multisimCommRec,
                 const MDLogger&   logger,
                 gmx_output_env_t* outputEnv) :
        fplog_{ fplog },
        commRec_{ commRec },
        multisimCommRec_{ multisimCommRec },
        logger_{ logger },
        outputEnv_{ outputEnv }
    {
    }

    //! Handle to log file.
    FILE* fplog_;
    //! Handle to communication record.
    t_commrec* commRec_;
    //! Handle to multisim communication record.
    const gmx_multisim_t* multisimCommRec_;
    //! Handle to propper logging framework.
    const MDLogger& logger_;
    //! Handle to file output handling.
    const gmx_output_env_t* outputEnv_;
};

/*! \brief
 * Collection of profiling information.
 */
class Profiling
{
public:
    //! Build profiling information collection.
    Profiling(t_nrnb* nrnb, gmx_walltime_accounting* walltimeAccounting, gmx_wallcycle* wallCycle) :
        nrnb(nrnb),
        wallCycle(wallCycle),
        walltimeAccounting(walltimeAccounting)
    {
    }

    //! Handle to datastructure.
    t_nrnb* nrnb;
    //! Handle to wallcycle stuff.
    gmx_wallcycle* wallCycle;
    //! Handle to wallcycle time accounting stuff.
    gmx_walltime_accounting* walltimeAccounting;
};

/*! \brief
 * Collection of constraint parameters.
 */
class ConstraintsParam
{
public:
    //! Build collection with handle to actual objects.
    ConstraintsParam(Constraints* constraints, gmx_enfrot* enforcedRotation, VirtualSitesHandler* vSite) :
        constr(constraints),
        enforcedRotation(enforcedRotation),
        vsite(vSite)
    {
    }

    //! Handle to constraint object.
    Constraints* constr;
    //! Handle to information about using enforced rotation.
    gmx_enfrot* enforcedRotation;
    //! Handle to vsite stuff.
    VirtualSitesHandler* vsite;
};

/*! \brief
 * Collection of legacy input information.
 */
class LegacyInput
{
public:
    //! Build collection from legacy input data.
    LegacyInput(int filenamesSize, const t_filenm* filenamesData, t_inputrec* inputRec, t_forcerec* forceRec) :
        numFile(filenamesSize),
        filenames(filenamesData),
        inputrec(inputRec),
        forceRec(forceRec)
    {
    }

    //! Number of input files.
    int numFile;
    //! File names.
    const t_filenm* filenames;
    //! Handle to simulation input record.
    t_inputrec* inputrec;
    //! Handle to simulation force record.
    t_forcerec* forceRec;
};

/*! \brief SimulatorBuilder parameter type for InteractiveMD.
 *
 * Conveys a non-owning pointer to implementation details.
 *
 * \todo If adding doxygen stubs actual add the full stub.
 */
class InteractiveMD
{
public:
    //! Create handle to IMD information.
    explicit InteractiveMD(ImdSession* imdSession) : imdSession(imdSession) {}

    //! Internal handle to IMD info.
    ImdSession* imdSession;
};

class SimulatorModules
{
public:
    SimulatorModules(IMDOutputProvider* mdOutputProvider, const MdModulesNotifier& notifier) :
        outputProvider(mdOutputProvider),
        mdModulesNotifier(notifier)
    {
    }

    IMDOutputProvider*       outputProvider;
    const MdModulesNotifier& mdModulesNotifier;
};

class CenterOfMassPulling
{
public:
    explicit CenterOfMassPulling(pull_t* pullWork) : pull_work(pullWork) {}

    pull_t* pull_work;
};

/*! \brief
 * Parameter type for IonSwapping SimulatorBuilder component.
 *
 * Conveys a non-owning pointer to implementation details.
 *
 * \todo Add full information.
 */
class IonSwapping
{
public:
    //! Create handle.
    IonSwapping(t_swap* ionSwap) : ionSwap(ionSwap) {}

    //! Internal storage for handle.
    t_swap* ionSwap;
};

/*! \brief
 * Collection of handles to topology information.
 */
class TopologyData
{
public:
    //! Build collection from simulation data.
    TopologyData(gmx_mtop_t* globalTopology, MDAtoms* mdAtoms) :
        top_global(globalTopology),
        mdAtoms(mdAtoms)
    {
    }

    //! Handle to global simulation topology.
    gmx_mtop_t* top_global;
    //! Handle to information about MDAtoms.
    MDAtoms* mdAtoms;
};

/*! \brief
 * Handle to information about the box.
 *
 * Design note: The client may own the BoxDeformation via std::unique_ptr, but we are not
 * transferring ownership at this time. (May be the subject of future changes.)
 */
class BoxDeformationHandle
{
public:
    //! Build handle to box stuff.
    BoxDeformationHandle(BoxDeformation* boxDeformation) : deform(boxDeformation) {}

    //! Internal storage for handle.
    BoxDeformation* deform;
};

/*! \libinternal
 * \brief Class preparing the creation of Simulator objects
 *
 * Objects of this class build Simulator objects, which in turn are used to
 * run molecular simulations.
 */
class SimulatorBuilder
{
public:
    void add(MembedHolder&& membedHolder);

    void add(std::unique_ptr<StopHandlerBuilder> stopHandlerBuilder)
    {
        stopHandlerBuilder_ = std::move(stopHandlerBuilder);
    }

    void add(SimulatorStateData&& simulatorStateData)
    {
        simulatorStateData_ = std::make_unique<SimulatorStateData>(simulatorStateData);
    }

    void add(SimulatorConfig&& simulatorConfig)
    {
        // Note: SimulatorConfig appears to the compiler to be trivially copyable,
        // but this may not be safe and may change in the future.
        simulatorConfig_ = std::make_unique<SimulatorConfig>(simulatorConfig);
    }

    void add(SimulatorEnv&& simulatorEnv)
    {
        simulatorEnv_ = std::make_unique<SimulatorEnv>(simulatorEnv);
    }

    void add(Profiling&& profiling) { profiling_ = std::make_unique<Profiling>(profiling); }

    void add(ConstraintsParam&& constraintsParam)
    {
        constraintsParam_ = std::make_unique<ConstraintsParam>(constraintsParam);
    }

    void add(LegacyInput&& legacyInput)
    {
        legacyInput_ = std::make_unique<LegacyInput>(legacyInput);
    }

    void add(ReplicaExchangeParameters&& replicaExchangeParameters);

    void add(InteractiveMD&& interactiveMd)
    {
        interactiveMD_ = std::make_unique<InteractiveMD>(interactiveMd);
    }

    void add(SimulatorModules&& simulatorModules)
    {
        simulatorModules_ = std::make_unique<SimulatorModules>(simulatorModules);
    }

    void add(CenterOfMassPulling&& centerOfMassPulling)
    {
        centerOfMassPulling_ = std::make_unique<CenterOfMassPulling>(centerOfMassPulling);
    }

    void add(IonSwapping&& ionSwapping)
    {
        ionSwapping_ = std::make_unique<IonSwapping>(ionSwapping);
    }

    void add(TopologyData&& topologyData)
    {
        topologyData_ = std::make_unique<TopologyData>(topologyData);
    }

    void add(BoxDeformationHandle&& boxDeformation)
    {
        boxDeformation_ = std::make_unique<BoxDeformationHandle>(boxDeformation);
    }

    /*!
     * \brief Pass the read checkpoint data for modular simulator
     *
     * Note that this is currently the point at which the ReadCheckpointDataHolder
     * is fully filled. Consequently it stops being an object at which read
     * operations from file are targeted, and becomes a read-only object from
     * which elements read their data to recreate an earlier internal state.
     *
     * Currently, this behavior change is not enforced. Once input reading and
     * simulator builder have matured, these restrictions could be imposed.
     *
     * See #3656
     */
    void add(std::unique_ptr<ReadCheckpointDataHolder> modularSimulatorCheckpointData);

    /*! \brief Build a Simulator object based on input data
     *
     * Return a pointer to a simulation object. The use of a parameter
     * pack insulates the builder from changes to the arguments of the
     * Simulator objects.
     *
     * \throws gmx::APIError if expected set-up methods have not been called before build()
     *
     * \return  Unique pointer to a Simulator object
     */
    std::unique_ptr<ISimulator> build(bool useModularSimulator);

private:
    // Note: we use std::unique_ptr instead of std::optional because we want to
    // allow for opaque types at the discretion of the module developer.
    /*! \brief Collection of handles to  individual information. */
    /*! \{ */
    std::unique_ptr<SimulatorConfig>           simulatorConfig_;
    std::unique_ptr<MembedHolder>              membedHolder_;
    std::unique_ptr<StopHandlerBuilder>        stopHandlerBuilder_;
    std::unique_ptr<SimulatorStateData>        simulatorStateData_;
    std::unique_ptr<SimulatorEnv>              simulatorEnv_;
    std::unique_ptr<Profiling>                 profiling_;
    std::unique_ptr<ConstraintsParam>          constraintsParam_;
    std::unique_ptr<LegacyInput>               legacyInput_;
    std::unique_ptr<ReplicaExchangeParameters> replicaExchangeParameters_;
    std::unique_ptr<InteractiveMD>             interactiveMD_;
    std::unique_ptr<SimulatorModules>          simulatorModules_;
    std::unique_ptr<CenterOfMassPulling>       centerOfMassPulling_;
    std::unique_ptr<IonSwapping>               ionSwapping_;
    std::unique_ptr<TopologyData>              topologyData_;
    std::unique_ptr<BoxDeformationHandle>      boxDeformation_;
    //! Contains checkpointing data for the modular simulator
    std::unique_ptr<ReadCheckpointDataHolder> modularSimulatorCheckpointData_;
    /*! \} */
};

} // namespace gmx

#endif // GMX_MDRUN_SIMULATORBUILDER_SIMULATORBUILDER_H
