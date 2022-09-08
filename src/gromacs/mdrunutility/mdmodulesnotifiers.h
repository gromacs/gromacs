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
/*! \libinternal \file
 * \brief
 * Declares gmx::MDModulesNotifiers.
 *
 * \author Christian Blau <blau@kth.se>
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */

#ifndef GMX_MDRUNUTILITY_MDMODULESNOTIFIERS_H
#define GMX_MDRUNUTILITY_MDMODULESNOTIFIERS_H

#include <string>
#include <vector>

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdrunutility/mdmodulesnotifier.h"

struct t_commrec;
struct gmx_mtop_t;
class WarningHandler;
enum class PbcType : int;

namespace gmx
{

class KeyValueTreeObject;
class KeyValueTreeObjectBuilder;
class LocalAtomSetManager;
class MDLogger;
class IndexGroupsAndNames;
class SeparatePmeRanksPermitted;
struct MDModulesCheckpointReadingDataOnMaster;
struct MDModulesCheckpointReadingBroadcast;
struct MDModulesWriteCheckpointData;

/*! \libinternal \brief Check if module outputs energy to a specific field.
 *
 * Ensures that energy is output for this module.
 */
struct MDModulesEnergyOutputToDensityFittingRequestChecker
{
    //! Trigger output to density fitting energy field
    bool energyOutputToDensityFitting_ = false;
};

/*! \libinternal \brief Check if QMMM module outputs energy to a specific field.
 *
 * Ensures that energy is output for QMMM module.
 */
struct MDModulesEnergyOutputToQMMMRequestChecker
{
    //! Trigger output to density fitting energy field
    bool energyOutputToQMMM_ = false;
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

/*! \libinternal \brief Provides coordinates and simulation box.
 */
struct CoordinatesAndBoxPreprocessed
{
    ArrayRefWithPadding<RVec> coordinates_;
    matrix                    box_;
    PbcType                   pbc_;
};

/*! \libinternal \brief Mdrun input filename.
 */
struct MdRunInputFilename
{
    //! The name of the run input file (.tpr) as output by grompp
    std::string mdRunFilename_;
};

/*! \libinternal \brief Notification for QM program input filename
 *  provided by user as command-line argument for grompp
 */
struct QMInputFileName
{
    //! Flag if QM Input File has been provided by user
    bool hasQMInputFileName_ = false;
    //! The name of the QM Input file (.inp)
    std::string qmInputFileName_;
};

/*! \libinternal
 * \brief Group of notifers to organize that MDModules
 * can receive callbacks they subscribe to.
 *
 * MDModules use members of this struct to subscribe to notifications
 * of particular events. When the event occurs, the callback provided
 * by a particular MDModule will be passed a parameter of the
 * particular type they are interested in.
 *
 * Typically, during the setup phase, modules subscribe to notifiers
 * that interest them by passing callbacks that expect a single parameter
 * that describes the event. These are stored for later use. See the
 * sequence diagram that follows:
   \msc
wordwraparcs=true,
hscale="2";

modules [label = "mdModules:\nMDModules"],
notifiers [label="notifiers\nMDModulesNotifiers"],
notifier [label="exampleNotifier:\nBuildMDModulesNotifier\n<EventX, EventY>::type"],
moduleA [label="moduleA"],
moduleB [label="moduleB"],
moduleC [label="moduleC"];

modules box moduleC [label = "mdModules creates and owns moduleA, moduleB, and moduleC"];
modules =>> notifiers [label="creates"];
notifiers =>> notifier [label="creates"];
notifier =>> notifiers [label="returns"];
notifiers =>> modules [label="returns"];

modules =>> moduleA [label="provides notifiers"];
moduleA =>> moduleA [label="unpacks\nnotifiers.exampleNotifier"];
moduleA =>> notifier [label="subscribes with\ncallback(EventX&)"];
notifier =>> notifier [label="records subscription\nto EventX"];
moduleA =>> notifier [label="subscribes with\ncallback(EventY&)"];
notifier =>> notifier [label="records subscription\nto EventY"];
moduleA =>> modules [label="returns"];

modules =>> moduleB [label="provides notifiers"];
moduleB =>> moduleB [label="unpacks\nnotifiers.exampleNotifier"];
moduleA =>> notifier [label="subscribes with\ncallback(EventY&)"];
notifier =>> notifier [label="records subscription\nto EventY"];
moduleB =>> modules [label="returns"];

modules =>> moduleC [label="provides notifiers"];
moduleC =>> moduleC [label="unpacks and keeps\nnotifiers.exampleNotifier"];
moduleC =>> modules [label="returns"];

   \endmsc

   * When the event occurs later on, the stored callbacks are used to
   * allow the modules to react. See the following sequence diagram,
   * which assumes that exampleNotifier was configured as in the
   * previous sequence diagram.

   \msc
wordwraparcs=true,
hscale="2";

moduleC [label="moduleC"],
notifier [label="exampleNotifier:\nBuildMDModulesNotifier\n<EventX, EventY>::type"],
moduleA [label="moduleA"],
moduleB [label="moduleB"];

moduleC box moduleB [label = "Later, when ModuleC is doing work"];
moduleC =>> moduleC [label="generates EventX"];
moduleC =>> moduleC [label="generates EventY"];
moduleC =>> notifier [label="calls notify(eventX)"];
notifier =>> moduleA [label="calls callback(eventX)"];
moduleA =>> moduleA [label="reacts to eventX"];
moduleA =>> notifier [label="returns"];

notifier =>> moduleC [label="returns"];
moduleC =>> notifier [label="calls notify(eventY)"];
notifier =>> moduleA [label="calls callback(eventY)"];
moduleA =>> moduleA [label="reacts to eventY"];
moduleA =>> notifier [label="returns"];
notifier =>> moduleB [label="calls callback(eventY)"];
moduleB =>> moduleB [label="reacts to eventY"];
moduleB =>> notifier [label="returns"];
notifier =>> moduleC [label="returns"];
   \endmsc
 *
 * The template arguments to the members of this struct are the
 * parameters passed to the callback functions, one type per
 * callback. Arguments passed as pointers are always meant to be
 * modified, but never meant to be stored (in line with the policy
 * everywhere else).
 *
 */
struct MDModulesNotifiers
{
    /*! \brief Pre-processing callback functions.
     * CoordinatesAndBoxPreprocessed Allows modules to access coordinates,
     *                                box and pbc during grompp
     * MDLogger Allows MdModule to use standard logging class for messages output
     * warninp* Allows modules to make grompp warnings, notes and errors
     * EnergyCalculationFrequencyErrors* allows modules to check if they match
     *                                   their required calculation frequency
     *                                   and add their error message if needed
     *                                   to the collected error messages
     * gmx_mtop_t* Allows modules to modify the topology during pre-processing
     * IndexGroupsAndNames provides modules with atom indices and their names
     * KeyValueTreeObjectBuilder enables writing of module internal data to
     *                           .tpr files.
     * QMInputFileName Allows QMMM module to know if user provided external QM input file
     */
    BuildMDModulesNotifier<const CoordinatesAndBoxPreprocessed&,
                           const MDLogger&,
                           WarningHandler*,
                           EnergyCalculationFrequencyErrors*,
                           gmx_mtop_t*,
                           const IndexGroupsAndNames&,
                           KeyValueTreeObjectBuilder,
                           const QMInputFileName&>::type preProcessingNotifier_;

    /*! \brief Handles subscribing and calling checkpointing callback functions.
     *
     * MDModulesCheckpointReadingDataOnMaster provides modules with their
     *                                        checkpointed data on the master
     *                                        node and checkpoint file version
     * MDModulesCheckpointReadingBroadcast provides modules with a communicator
     *                                     and the checkpoint file version to
     *                                     distribute their data
     * MDModulesWriteCheckpointData provides the modules with a key-value-tree
     *                              builder to store their checkpoint data and
     *                              the checkpoint file version
     */
    BuildMDModulesNotifier<MDModulesCheckpointReadingDataOnMaster, MDModulesCheckpointReadingBroadcast, MDModulesWriteCheckpointData>::type
            checkpointingNotifier_;

    /*! \brief Handles subscribing and calling callbacks during simulation setup.
     *
     * const KeyValueTreeObject& provides modules with the internal data they
     *                           wrote to .tpr files
     * LocalAtomSetManager* enables modules to add atom indices to local atom sets
     *                      to be managed
     * const MDLogger& Allows MdModule to use standard logging class for messages output
     * const gmx_mtop_t& provides the topology of the system to the modules
     * MDModulesEnergyOutputToDensityFittingRequestChecker* enables modules to
     *                      report if they want to write their energy output
     *                      to the density fitting field in the energy files
     * MDModulesEnergyOutputToQMMMRequestChecker* enables QMMM module to
     *                      report if it want to write their energy output
     *                      to the "Quantum En." field in the energy files
     * SeparatePmeRanksPermitted* enables modules to report if they want
     *                      to disable dedicated PME ranks
     * const PbcType& provides modules with the periodic boundary condition type
     *                that is used during the simulation
     * const SimulationTimeStep& provides modules with the simulation time-step
     *                           that allows them to interconvert between step
     *                           time information
     * const t_commrec& provides a communicator to the modules during simulation
     *                  setup
     * const MdRunInputFilename& Allows modules to know .tpr filename during mdrun
     */
    BuildMDModulesNotifier<const KeyValueTreeObject&,
                           LocalAtomSetManager*,
                           const MDLogger&,
                           const gmx_mtop_t&,
                           MDModulesEnergyOutputToDensityFittingRequestChecker*,
                           MDModulesEnergyOutputToQMMMRequestChecker*,
                           SeparatePmeRanksPermitted*,
                           const PbcType&,
                           const SimulationTimeStep&,
                           const t_commrec&,
                           const MdRunInputFilename&>::type simulationSetupNotifier_;
};

} // namespace gmx

#endif
