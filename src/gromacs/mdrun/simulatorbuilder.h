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

#include "gromacs/mdlib/vsite.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/mdmodulenotification.h"

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
class ImdSession;
class MDLogger;
class MDAtoms;
class ISimulator;
class StopHandlerBuilder;
struct MdrunOptions;

/*! \brief Membed SimulatorBuilder parameter type.
 *
 * Does not (yet) encapsulate ownership semantics of resources. Simulator is
 * not (necessarily) granted ownership of resources. Client is responsible for
 * maintaining the validity of resources for the life time of the Simulator,
 * then for cleaning up those resources.
 */
class MembedHolder
{
public:
    explicit MembedHolder(gmx_membed_t* membed) : membed_(membed) {}

    gmx_membed_t* membed() { return membed_; }

private:
    gmx_membed_t* membed_;
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

/*! \libinternal
 * \brief Class preparing the creation of Simulator objects
 *
 * Objects of this class build Simulator objects, which in turn are used to
 * run molecular simulations.
 */
class SimulatorBuilder
{
public:
    void add(MembedHolder&& membedHolder)
    {
        membedHolder_ = std::make_unique<MembedHolder>(membedHolder);
    }

    void add(std::unique_ptr<StopHandlerBuilder> stopHandlerBuilder)
    {
        stopHandlerBuilder_ = std::move(stopHandlerBuilder);
    }

    void add(SimulatorStateData&& simulatorStateData)
    {
        simulatorStateData_ = std::make_unique<SimulatorStateData>(simulatorStateData);
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
    std::unique_ptr<ISimulator> build(bool                             useModularSimulator,
                                      FILE*                            fplog,
                                      t_commrec*                       cr,
                                      const gmx_multisim_t*            ms,
                                      const MDLogger&                  mdlog,
                                      int                              nfile,
                                      const t_filenm*                  fnm,
                                      const gmx_output_env_t*          oenv,
                                      const MdrunOptions&              mdrunOptions,
                                      StartingBehavior                 startingBehavior,
                                      VirtualSitesHandler*             vsite,
                                      Constraints*                     constr,
                                      gmx_enfrot*                      enforcedRotation,
                                      BoxDeformation*                  deform,
                                      IMDOutputProvider*               outputProvider,
                                      const MdModulesNotifier&         mdModulesNotifier,
                                      t_inputrec*                      inputrec,
                                      ImdSession*                      imdSession,
                                      pull_t*                          pull_work,
                                      t_swap*                          swap,
                                      gmx_mtop_t*                      top_global,
                                      MDAtoms*                         mdAtoms,
                                      t_nrnb*                          nrnb,
                                      gmx_wallcycle*                   wcycle,
                                      t_forcerec*                      fr,
                                      MdrunScheduleWorkload*           runScheduleWork,
                                      const ReplicaExchangeParameters& replExParams,
                                      gmx_walltime_accounting*         walltime_accounting,
                                      bool                             doRerun);

private:
    std::unique_ptr<MembedHolder>       membedHolder_;
    std::unique_ptr<StopHandlerBuilder> stopHandlerBuilder_;
    std::unique_ptr<SimulatorStateData> simulatorStateData_;
};

} // namespace gmx

#endif // GMX_MDRUN_SIMULATORBUILDER_SIMULATORBUILDER_H
