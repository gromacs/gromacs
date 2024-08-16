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

#include <cstdint>

#include <optional>
#include <string>
#include <vector>

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/matrix.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdrunutility/mdmodulesnotifier.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"


struct t_commrec;
struct gmx_mtop_t;
class WarningHandler;
enum class PbcType : int;
struct t_inputrec;

namespace gmx
{

class KeyValueTreeObject;
class KeyValueTreeObjectBuilder;
class LocalAtomSetManager;
class MDLogger;
class IndexGroupsAndNames;
class SeparatePmeRanksPermitted;
struct MDModulesCheckpointReadingDataOnMain;
struct MDModulesCheckpointReadingBroadcast;
struct MDModulesWriteCheckpointData;
enum class StartingBehavior;

/*! \libinternal \brief Notification that atoms may have been redistributed
 *
 * This notification is emitted at the end of the DD (re)partioning
 * or without DD right after atoms have put into the box.
 * The local atom sets are updated for the new atom order when this signal is emitted.
 * The coordinates of atoms can be shifted by periodic vectors
 * before the signal was emitted.
 */
struct MDModulesAtomsRedistributedSignal
{
    MDModulesAtomsRedistributedSignal(const matrix box, gmx::ArrayRef<const RVec> x) :
        box_(createMatrix3x3FromLegacyMatrix(box)), x_(x)
    {
    }

    //! The simulation unit cell
    const Matrix3x3 box_;
    //! List of local atom coordinates after partitioning
    gmx::ArrayRef<const RVec> x_;
};

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
    EnergyCalculationFrequencyErrors(int64_t energyCalculationIntervalInSteps);
    //! Return the number of steps of an energy calculation interval
    std::int64_t energyCalculationIntervalInSteps() const;
    //! Collect error messages
    void addError(const std::string& errorMessage);
    //! Return error messages
    const std::vector<std::string>& errorMessages() const;

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

/*! \libinternal \brief Energy trajectory output filename from Mdrun.
 */
struct EdrOutputFilename
{
    //! The name of energy output file
    std::string edrOutputFilename_;
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

/*! \libinternal \brief Provides the constant ensemble temperature
 */
struct EnsembleTemperature
{
    /*! \libinternal
     * \brief Check whether the constant ensemble temperature is available.
     * Then, store the value as optional.
     */
    explicit EnsembleTemperature(const t_inputrec& ir);
    //! The constant ensemble temperature
    std::optional<real> constantEnsembleTemperature_;
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
     *
     * \tparam CoordinatesAndBoxPreprocessed
     *                              Allows modules to access coordinates,
     *                              box and PBC during grompp
     * \tparam MDLogger             Allows MdModule to use standard logging class for
     *                              output of messages
     * \tparam warninp*             Allows modules to make grompp warnings, notes and errors
     * \tparam EnergyCalculationFrequencyErrors*
     *                              Allows modules to check if they match their required calculation
     *                              frequency and add their error message if needed to the
     *                              collected error messages
     * \tparam gmx_mtop_t*          Allows modules to modify the topology during pre-processing
     * \tparam IndexGroupsAndNames  Provides modules with atom indices and their names
     * \tparam KeyValueTreeObjectBuilder
     *                              Enables writing of module internal data to .tpr files.
     * \tparam QMInputFileName      Allows the QMMM module to know if the user has provided
     *                              an external QM input file
     * \tparam EnsembleTemperature  Provides modules with the constant ensemble temperature.
     */
    BuildMDModulesNotifier<const CoordinatesAndBoxPreprocessed&,
                           const MDLogger&,
                           WarningHandler*,
                           EnergyCalculationFrequencyErrors*,
                           gmx_mtop_t*,
                           const IndexGroupsAndNames&,
                           KeyValueTreeObjectBuilder,
                           const QMInputFileName&,
                           const EnsembleTemperature&>::type preProcessingNotifier_;

    /*! \brief Handles subscribing and calling checkpointing callback functions.
     *
     * \tparam MDModulesCheckpointReadingDataOnMain
     *                              Provides modules with their checkpointed data
     *                              on the main node and checkpoint file version
     * \tparam MDModulesCheckpointReadingBroadcast
     *                              Provides modules with a communicator and the
     *                              checkpoint file version to distribute their data
     * \tparam MDModulesWriteCheckpointData
     *                              Provides the modules with a key-value-tree
     *                              builder to store their checkpoint data and
     *                              the checkpoint file version
     */
    BuildMDModulesNotifier<MDModulesCheckpointReadingDataOnMain, MDModulesCheckpointReadingBroadcast, MDModulesWriteCheckpointData>::type
            checkpointingNotifier_;

    /*! \brief Handles subscribing and calling callbacks during simulation setup.
     *
     * \tparam KeyValueTreeObject&  Provides modules with the internal data they
     *                              wrote to .tpr files
     * \tparam LocalAtomSetManager* Enables modules to add atom indices to local atom sets
     *                              to be managed
     * \tparam StartingBehavio&     Provides modules with the starting behavior of the simulation
     * \tparam MDLogger&            Allows MdModule to use standard logging class for messages
     *                              output
     * \tparam gmx_mtop_t&          Provides the topology of the system to the modules
     * \tparam MDModulesEnergyOutputToDensityFittingRequestChecker*
     *                              Enables modules to report if they want to write their
     *                              energy output to the density fitting field in the energy files
     * \tparam MDModulesEnergyOutputToQMMMRequestChecker*
     *                              Enables QMMM module to report if it wants to write its energy
     *                              output to the "Quantum En." field in the energy files
     * \tparam SeparatePmeRanksPermitted*
     *                              Enables modules to report if they want to disable dedicated
     *                              PME ranks
     * \tparam PbcType&             Provides modules with the periodic boundary condition type
     *                              that is used during the simulation
     * \tparam SimulationTimeStep&  Provides modules with the simulation time-step that allows
     *                              them to interconvert between step and time information
     * \tparam EnsembleTemperature& Provides modules with the (eventual) constant ensemble
     *                              temperature
     * \tparam t_commrec&           Provides a communicator to the modules during simulation
     *                              setup
     * \tparam MdRunInputFilename&  Allows modules to know .tpr filename during mdrun
     * \tparam EdrOutputFilename&   Allows modules to know .edr filename during mdrun
     */
    BuildMDModulesNotifier<const KeyValueTreeObject&,
                           LocalAtomSetManager*,
                           const StartingBehavior&,
                           const MDLogger&,
                           const gmx_mtop_t&,
                           const MDModulesAtomsRedistributedSignal,
                           MDModulesEnergyOutputToDensityFittingRequestChecker*,
                           MDModulesEnergyOutputToQMMMRequestChecker*,
                           SeparatePmeRanksPermitted*,
                           const PbcType&,
                           const SimulationTimeStep&,
                           const EnsembleTemperature&,
                           const t_commrec&,
                           const MdRunInputFilename&,
                           const EdrOutputFilename&>::type simulationSetupNotifier_;
};

} // namespace gmx

#endif
