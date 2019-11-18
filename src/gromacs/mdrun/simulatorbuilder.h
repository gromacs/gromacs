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
/*! \internal
 * \brief Declares the simulator builder for mdrun
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_SIMULATORBUILDER_H
#define GMX_MDRUN_SIMULATORBUILDER_H

#include <memory>

#include "gromacs/modularsimulator/modularsimulator.h"

#include "legacysimulator.h"

class energyhistory_t;
struct gmx_enerdata_t;
struct gmx_enfrot;
struct gmx_mtop_t;
struct gmx_membed_t;
struct gmx_multisim_t;
struct gmx_output_env_t;
struct gmx_vsite_t;
struct gmx_wallcycle;
struct gmx_walltime_accounting;
struct ObservablesHistory;
struct pull_t;
struct ReplicaExchangeParameters;
struct t_commrec;
struct t_fcdata;
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

/*! \libinternal
 * \brief Class preparing the creation of Simulator objects
 *
 * Objects of this class build Simulator objects, which in turn are used to
 * run molecular simulations. Currently, this only has a single public
 * `build` function which takes all arguments needed to build the
 * `LegacySimulator`.
 */
class SimulatorBuilder
{
public:
    /*! \brief Build a Simulator object based on input data
     *
     * Return a pointer to a simulation object. The use of a parameter
     * pack insulates the builder from changes to the arguments of the
     * Simulator objects.
     *
     * @return  Unique pointer to a Simulator object
     */
    template<typename... Args>
    std::unique_ptr<ISimulator> build(bool inputIsCompatibleWithModularSimulator, Args&&... args);
};


//! Build a Simulator object
template<typename... Args>
std::unique_ptr<ISimulator> SimulatorBuilder::build(bool inputIsCompatibleWithModularSimulator, Args&&... args)
{
    // GMX_DISABLE_MODULAR_SIMULATOR allows to disable modular simulator for all uses
    const auto disableModularSimulator = (getenv("GMX_DISABLE_MODULAR_SIMULATOR") != nullptr);

    if (!disableModularSimulator && inputIsCompatibleWithModularSimulator)
    {
        // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
        return std::unique_ptr<ModularSimulator>(new ModularSimulator(std::forward<Args>(args)...));
    }
    // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
    return std::unique_ptr<LegacySimulator>(new LegacySimulator(std::forward<Args>(args)...));
}

} // namespace gmx

#endif // GMX_MDRUN_SIMULATORBUILDER_SIMULATORBUILDER_H
