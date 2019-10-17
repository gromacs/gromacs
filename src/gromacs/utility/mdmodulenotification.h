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
/*! \libinternal \file
 * \brief
 * Declares gmx::MdModuleNotification.
 *
 * \author Christian Blau <blau@kth.se>
 * \inlibraryapi
 * \ingroup module_utility
 */

#ifndef GMX_MDRUNUTILITY_MDMODULENOTIFICATION_H
#define GMX_MDRUNUTILITY_MDMODULENOTIFICATION_H

#include <functional>
#include <string>
#include <vector>

struct t_commrec;

namespace gmx
{

/*! \libinternal \brief
 * Subscribe and trigger notification functions.
 *
 * Extends MdModuleNotificationBase with new notification function and routine
 * to subscribe new listeners.
 *
 * To create a class of this type that provides callbacks, e.g., for events
 * EventA, and EventB use registerMdModuleNotification<EventA, EventB>::type.
 *
 * \tparam CallParameter of the function to be notified
 * \tparam MdModuleNotificationBase class to be extended with a notification
 *                                  with CallParameter
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
 * \note All added subscribers are required to out-live the MdModuleNotification
 *
 */
template<class CallParameter, class MdModuleNotificationBase>
class MdModuleNotification : public MdModuleNotificationBase
{
public:
    //! Make base class notification trigger available to this class
    using MdModuleNotificationBase::notify;
    //! Make base class subscription available to this class
    using MdModuleNotificationBase::subscribe;

    /*! \brief Trigger the subscribed notifications.
     * \param[in] callParameter of the function to be called back
     */
    void notify(CallParameter callParameter) const
    {
        for (auto& callBack : callBackFunctions_)
        {
            callBack(callParameter);
        }
    }

    /*! \brief
     * Add callback function to be called when notification is triggered.
     *
     * Notifications are distinguished by their call signature.
     *
     * \param[in] callBackFunction to be called from this class
     */
    void subscribe(std::function<void(CallParameter)> callBackFunction)
    {
        callBackFunctions_.emplace_back(callBackFunction);
    }

private:
    std::vector<std::function<void(CallParameter)>> callBackFunctions_;
};

/*! \internal
 * \brief Aide to avoid nested MdModuleNotification definition.
 *
 * Instead of
 * MdModuleNotification<CallParameterA, MdModuleNotification<CallParameterB, etc ... >>
 * this allows to write
 * registerMdModuleNotification<CallParameterA, CallParameterB, ...>::type
 *
 * \tparam CallParameter all the event types to be registered
 */
template<class... CallParameter>
struct registerMdModuleNotification;

/*! \internal \brief Template specialization to end parameter unpacking recursion.
 */
template<>
struct registerMdModuleNotification<>
{
    /*! \internal
     * \brief Do nothing but be base class of MdModuleNotification.
     *
     * Required so that using MdModuleNotificationBase::notify and
     * MdModuleNotificationBase::subscribe are valid in derived class.
     */
    class NoCallParameter
    {
    public:
        //! Do nothing but provide MdModuleNotification::notify to derived class
        void notify() {}
        //! Do nothing but provide MdModuleNotification::subscribe to derived class
        void subscribe() {}
    };
    /*! \brief Defines a type if no notifications are managed.
     *
     * This ensures that code works with MdModuleCallParameterManagement that
     * does not manage any notifications.
     */
    using type = NoCallParameter;
};

/*! \libinternal
 * \brief Template specialization to assemble MdModuleNotification.
 *
 * Assembly of MdModuleNotification is performed by recursively taking off the
 * front of the CallParameter parameter pack and constructing the nested type
 * definition of MdModuleNotification base classes.
 *
 * \tparam CurrentCallParameter front of the template parameter pack
 * \tparam CallParameter rest of the event types
 */
template<class CurrentCallParameter, class... CallParameter>
struct registerMdModuleNotification<CurrentCallParameter, CallParameter...>
{
    // private:
    //! The next type with rest of the arguments with the front parameter removed.
    using next_type = typename registerMdModuleNotification<CallParameter...>::type;
    //! The type of the MdModuleNotification
    using type = MdModuleNotification<CurrentCallParameter, next_type>;
};

class KeyValueTreeObject;
class KeyValueTreeObjectBuilder;
class LocalAtomSetManager;
class IndexGroupsAndNames;
struct MdModulesCheckpointReadingDataOnMaster;
struct MdModulesCheckpointReadingBroadcast;
struct MdModulesWriteCheckpointData;
struct PeriodicBoundaryConditionType
{
    int pbcType;
};

struct MdModulesEnergyOutputToDensityFittingRequestChecker
{
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

struct SimulationTimeStep
{
    //! Time step (ps)
    const double delta_t;
};

struct MdModulesNotifier
{
    //! Register callback function types for MdModule
    registerMdModuleNotification<const t_commrec&,
                                 EnergyCalculationFrequencyErrors*,
                                 IndexGroupsAndNames,
                                 KeyValueTreeObjectBuilder,
                                 const KeyValueTreeObject&,
                                 LocalAtomSetManager*,
                                 MdModulesEnergyOutputToDensityFittingRequestChecker*,
                                 MdModulesCheckpointReadingDataOnMaster,
                                 MdModulesCheckpointReadingBroadcast,
                                 MdModulesWriteCheckpointData,
                                 PeriodicBoundaryConditionType,
                                 const SimulationTimeStep&>::type notifier_;
};

} // namespace gmx

#endif
