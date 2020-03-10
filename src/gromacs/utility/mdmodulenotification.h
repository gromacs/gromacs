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
 * \brief
 * Declares gmx::MdModulesNotifier.
 *
 * \author Christian Blau <blau@kth.se>
 * \inlibraryapi
 * \ingroup module_utility
 */

#ifndef GMX_UTILITY_MDMODULENOTIFICATION_H
#define GMX_UTILITY_MDMODULENOTIFICATION_H

#include <string>
#include <vector>

#include "gromacs/utility/mdmodulenotification-impl.h"

struct t_commrec;
enum class PbcType : int;

namespace gmx
{

class KeyValueTreeObject;
class KeyValueTreeObjectBuilder;
class LocalAtomSetManager;
class IndexGroupsAndNames;
struct MdModulesCheckpointReadingDataOnMaster;
struct MdModulesCheckpointReadingBroadcast;
struct MdModulesWriteCheckpointData;

/*! \libinternal \brief Check if module outputs energy to a specific field.
 *
 * Ensures that energy is output for this module.
 */
struct MdModulesEnergyOutputToDensityFittingRequestChecker
{
    //! Trigger output to density fitting energy field
    bool energyOutputToDensityFitting_ = false;
};

/*! \libinternal
 * \brief Collect errors for the energy calculation frequency.
 *
 * Collect errors regarding energy calculation frequencies as strings that then
 * may be used to issue errors.
 *
 * \note The mdp option "nstcalcenergy" is altered after reading the .mdp input
 *       and only used in certain integrators, thus this class is to be used
 *       only after all these operations are done.
 */
class EnergyCalculationFrequencyErrors
{
public:
    //! Construct by setting the energy calculation frequency
    EnergyCalculationFrequencyErrors(int64_t energyCalculationIntervalInSteps) :
        energyCalculationIntervalInSteps_(energyCalculationIntervalInSteps)
    {
    }
    //! Return the number of steps of an energy calculation interval
    std::int64_t energyCalculationIntervalInSteps() const
    {
        return energyCalculationIntervalInSteps_;
    }
    //! Collect error messages
    void addError(const std::string& errorMessage) { errorMessages_.push_back(errorMessage); }
    //! Return error messages
    const std::vector<std::string>& errorMessages() const { return errorMessages_; }

private:
    //! The frequency of energy calculations
    const std::int64_t energyCalculationIntervalInSteps_;
    //! The error messages
    std::vector<std::string> errorMessages_;
};

/*! \libinternal \brief Provides the simulation time step in ps.
 */
struct SimulationTimeStep
{
    //! Time step (ps)
    const double delta_t;
};

/*! \libinternal
 * \brief Collection of callbacks to MDModules at differnt run-times.
 *
 * MDModules use members of this struct to sign up for callback functionality.
 *
 * The members of the struct represent callbacks at these run-times:
 *
 *  When pre-processing the simulation data
 *  When reading and writing check-pointing data
 *  When setting up simulation after reading in the tpr file
 *
   \msc
   wordwraparcs=true,
   hscale="2";

   runner [label="runner:\nMdrunner"],
   CallParameter [label = "eventA:\nCallParameter"],
   MOD [label = "mdModules_:\nMdModules"],
   ModuleA [label="moduleA"],
   ModuleB [label="moduleB"],
   MdModuleNotification [label="notifier_:\nMdModuleNotification"];

   MOD box MdModuleNotification [label = "mdModules_ owns notifier_ and moduleA/B"];
   MOD =>> ModuleA [label="instantiates(notifier_)"];
   ModuleA =>> MdModuleNotification [label="subscribe(otherfunc)"];
   ModuleA =>> MOD;
   MOD =>> ModuleB [label="instantiates(notifier_)"];
   ModuleB =>> MdModuleNotification [label="subscribe(func)"];
   ModuleB =>> MOD;
   runner =>> CallParameter [label="instantiate"];
   CallParameter =>> runner ;
   runner =>> MOD [label="notify(eventA)"];
   MOD =>> MdModuleNotification [label="notify(eventA)"];
   MdModuleNotification =>> ModuleA [label="notify(eventA)"];
   ModuleA -> ModuleA [label="func(eventA)"];
   MdModuleNotification =>> ModuleB [label="notify(eventA)"];
   ModuleB -> ModuleB [label="otherfunc(eventA)"];

   \endmsc
 *
 * The template arguments to the members of this struct directly reflect
 * the callback function signature. Arguments passed as pointers are always
 * meant to be modified, but never meant to be stored (in line with the policy
 * everywhere else).
 */
struct MdModulesNotifier
{
    /*! \brief Pre-processing callback functions.
     *
     * EnergyCalculationFrequencyErrors* allows modules to check if they match
     *                                   their required calculation frequency
     *                                   and add their error message if needed
     *                                   to the collected error messages
     * IndexGroupsAndNames provides modules with atom indices and their names
     * KeyValueTreeObjectBuilder enables writing of module internal data to
     *                           .tpr files.
     */
    registerMdModuleNotification<EnergyCalculationFrequencyErrors*, IndexGroupsAndNames, KeyValueTreeObjectBuilder>::type preProcessingNotifications_;

    /*! \brief Checkpointing callback functions.
     *
     * MdModulesCheckpointReadingDataOnMaster provides modules with their
     *                                        checkpointed data on the master
     *                                        node and checkpoint file version
     * MdModulesCheckpointReadingBroadcast provides modules with a communicator
     *                                     and the checkpoint file version to
     *                                     distribute their data
     * MdModulesWriteCheckpointData provides the modules with a key-value-tree
     *                              builder to store their checkpoint data and
     *                              the checkpoint file version
     */
    registerMdModuleNotification<MdModulesCheckpointReadingDataOnMaster,
                                 MdModulesCheckpointReadingBroadcast,
                                 MdModulesWriteCheckpointData>::type checkpointingNotifications_;

    /*! \brief Callbacks during simulation setup.
     *
     * const KeyValueTreeObject& provides modules with the internal data they
     *                           wrote to .tpr files
     * LocalAtomSetManager* enables modules to add atom indices to local atom sets
     *                      to be managed
     * MdModulesEnergyOutputToDensityFittingRequestChecker* enables modules to
     *                      report if they want to write their energy output
     *                      to the density fitting field in the energy files
     * const PbcType& provides modules with the periodic boundary condition type
     *                that is used during the simulation
     * const SimulationTimeStep& provides modules with the simulation time-step
     *                           that allows them to interconvert between step
     *                           time information
     * const t_commrec& provides a communicator to the modules during simulation
     *                  setup
     */
    registerMdModuleNotification<const KeyValueTreeObject&,
                                 LocalAtomSetManager*,
                                 MdModulesEnergyOutputToDensityFittingRequestChecker*,
                                 const PbcType&,
                                 const SimulationTimeStep&,
                                 const t_commrec&>::type simulationSetupNotifications_;
};

} // namespace gmx

#endif
