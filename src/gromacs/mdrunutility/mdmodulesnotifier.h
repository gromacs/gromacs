/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Declares gmx::MDModulesNotifier and builder.
 *
 * \author Christian Blau <blau@kth.se>
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */

#ifndef GMX_MDRUNUTILITY_MDMODULESNOTIFIER_H
#define GMX_MDRUNUTILITY_MDMODULESNOTIFIER_H

#include <functional>
#include <vector>

namespace gmx
{

/*! \libinternal
 * \brief Organizes notifications about an event of interest to modules.
 *
 * An object of this type permits modules to subscribe to the
 * corresponding event. The template types of this type encode what
 * information is available when the event occurs. Modules \c
 * subscribe() by providing a callback function that accepts a single
 * parameter of such an event type. The code that handles that event
 * has the responsibilty to call \c notify() afterwards. The
 * subscribed modules then receive the callback with the requested
 * event type as an argument.
 *
 * See gmx::MDModulesNotifiers for sequence diagrams for an example.
 *
 * This suits scenarios where several objects are built (or re-built)
 * and one or more modules need to know when one or more of such
 * objects are available (or updated), so they can adapt their
 * internal state accordingly. Examples include responding to loading
 * input data, or to changes related to a recurring process like
 * checkpointing or partitioning. The coupling between these modules
 * is now expressed indirectly. This improves the modularity and
 * testability of those modules.
 *
 * The implementation provides the necessary flexibility to be
 * parameterized with multiple event types and provide \c callback()
 * and \b notify() methods corresponding to each related event. This
 * is done by inheriting from a series of base classes, each of which
 * handles a single type of event. BuildMDModulesNotifier implements
 * the details. To create a class of this type that provides two
 * events with callbacks that receive respectively types TypeA and
 * TypeB, use BuildMDModulesNotifier<TypeA, TypeB>::type.
 *
 * \tparam CallParameter of the function to be notified
 * \tparam MDModulesNotifierBase class to be extended with a notification
 *                                  with CallParameter
 *
 * \note All added subscribers are required to out-live the MDModulesNotifier
 *
 */
template<class CallParameter, class MDModulesNotifierBase>
class MDModulesNotifier : public MDModulesNotifierBase
{
public:
    //! Make base class notification trigger available to this class
    using MDModulesNotifierBase::notify;
    //! Make base class subscription available to this class
    using MDModulesNotifierBase::subscribe;

    /*! \brief Notifies subscribers of the event described by \c
     * callbackParameter.
     *
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
     * Add callback function to be called when \c notify() is called
     *
     * \param[in] callBackFunction to be called
     */
    void subscribe(std::function<void(CallParameter)> callBackFunction)
    {
        callBackFunctions_.emplace_back(callBackFunction);
    }

private:
    std::vector<std::function<void(CallParameter)>> callBackFunctions_;
};

/*! \internal
 * \brief Aide to avoid nested MDModulesNotifier definition.
 *
 * Instead of
 * MDModulesNotifier<CallParameterA, MDModulesNotifier<CallParameterB, etc ... >>
 * this allows to write
 * BuildMDModulesNotifier<CallParameterA, CallParameterB, ...>::type
 *
 * \tparam CallParameter all the callback types to be registered
 */
template<class... CallParameter>
struct BuildMDModulesNotifier;

/*! \internal \brief Template specialization to end parameter unpacking recursion.
 */
template<>
struct BuildMDModulesNotifier<>
{
    /*! \internal
     * \brief Do nothing but be base class of MDModulesNotifier.
     *
     * Required so that using MDModulesNotifierBase::notify and
     * MDModulesNotifierBase::subscribe are valid in derived class.
     */
    class NoCallParameter
    {
    public:
        //! Do nothing but provide MDModulesNotifier::notify to derived class
        void notify() {}
        //! Do nothing but provide MDModulesNotifier::subscribe to derived class
        void subscribe() {}
    };
    /*! \brief Defines a type if no notifications are managed.
     *
     * This ensures that code works with MDModuleCallParameterManagement that
     * does not manage any notifications.
     */
    using type = NoCallParameter;
};

/*! \libinternal
 * \brief Template specialization to assemble MDModulesNotifier.
 *
 * Assembly of MDModulesNotifier is performed by recursively taking off the
 * front of the CallParameter parameter pack and constructing the nested type
 * definition of MDModulesNotifier base classes.
 *
 * \tparam CurrentCallParameter front of the template parameter pack
 * \tparam CallParameter rest of the callback types
 */
template<class CurrentCallParameter, class... CallParameter>
struct BuildMDModulesNotifier<CurrentCallParameter, CallParameter...>
{
    // private:
    //! The next type with rest of the arguments with the front parameter removed.
    using next_type = typename BuildMDModulesNotifier<CallParameter...>::type;
    //! The type of the MDModulesNotifier
    using type = MDModulesNotifier<CurrentCallParameter, next_type>;
};

} // namespace gmx

#endif
