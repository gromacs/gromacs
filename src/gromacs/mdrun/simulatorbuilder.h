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
enum class StartingBehavior;
class StopHandlerBuilder;
class VirtualSitesHandler;

struct SimulatorConfig
{
public:
    SimulatorConfig(const MdrunOptions&    mdrunOptions,
                    StartingBehavior       startingBehavior,
                    MdrunScheduleWorkload* runScheduleWork) :
        mdrunOptions_(mdrunOptions),
        startingBehavior_(startingBehavior),
        runScheduleWork_(runScheduleWork)
    {
    }
    // TODO: Specify copy and move semantics.

    const MdrunOptions&    mdrunOptions_;
    StartingBehavior       startingBehavior_;
    MdrunScheduleWorkload* runScheduleWork_;
};


// TODO: Reconsider the name.
struct SimulatorStateData
{
    t_state*            globalState_p;
    ObservablesHistory* observablesHistory_p;
    gmx_enerdata_t*     enerdata_p;
    gmx_ekindata_t*     ekindata_p;

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

    SimulatorStateData(const SimulatorStateData& simulatorStateData) = default;
};

class SimulatorEnv
{
public:
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

    FILE*                   fplog_;
    t_commrec*              commRec_;
    const gmx_multisim_t*   multisimCommRec_;
    const MDLogger&         logger_;
    const gmx_output_env_t* outputEnv_;
};

class Profiling
{
public:
    Profiling(t_nrnb* nrnb, gmx_walltime_accounting* walltimeAccounting, gmx_wallcycle* wallCycle) :
        nrnb(nrnb),
        wallCycle(wallCycle),
        walltimeAccounting(walltimeAccounting)
    {
    }

    t_nrnb*                  nrnb;
    gmx_wallcycle*           wallCycle;
    gmx_walltime_accounting* walltimeAccounting;
};

class ConstraintsParam
{
public:
    ConstraintsParam(Constraints* constraints, gmx_enfrot* enforcedRotation, VirtualSitesHandler* vSite) :
        constr(constraints),
        enforcedRotation(enforcedRotation),
        vsite(vSite)
    {
    }

    Constraints*         constr;
    gmx_enfrot*          enforcedRotation;
    VirtualSitesHandler* vsite;
};

class LegacyInput
{
public:
    LegacyInput(int filenamesSize, const t_filenm* filenamesData, t_inputrec* inputRec, t_forcerec* forceRec) :
        numFile(filenamesSize),
        filenames(filenamesData),
        inputrec(inputRec),
        forceRec(forceRec)
    {
    }
    int             numFile;
    const t_filenm* filenames;
    t_inputrec*     inputrec;
    t_forcerec*     forceRec;
};

/*! \brief SimulatorBuilder parameter type for InteractiveMD.
 *
 * Conveys a non-owning pointer to implementation details.
 */
class InteractiveMD
{
public:
    explicit InteractiveMD(ImdSession* imdSession) : imdSession(imdSession) {}

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

/*!
 * \brief Parameter type for IonSwapping SimulatorBuilder component.
 *
 * Conveys a non-owning pointer to implementation details.
 */
class IonSwapping
{
public:
    IonSwapping(t_swap* ionSwap) : ionSwap(ionSwap) {}
    t_swap* ionSwap;
};

class TopologyData
{
public:
    TopologyData(gmx_mtop_t* globalTopology, MDAtoms* mdAtoms) :
        top_global(globalTopology),
        mdAtoms(mdAtoms)
    {
    }
    gmx_mtop_t* top_global;
    MDAtoms*    mdAtoms;
};

// Design note: The client may own the BoxDeformation via std::unique_ptr, but we are not
// transferring ownership at this time. (Maybe be the subject of future changes.)
class BoxDeformationHandle
{
public:
    BoxDeformationHandle(BoxDeformation* boxDeformation) : deform(boxDeformation) {}
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
};

} // namespace gmx

#endif // GMX_MDRUN_SIMULATORBUILDER_SIMULATORBUILDER_H
