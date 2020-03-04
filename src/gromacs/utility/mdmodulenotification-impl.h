/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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

#ifndef GMX_UTILITY_MDMODULENOTIFICATION_IMPL_H
#define GMX_UTILITY_MDMODULENOTIFICATION_IMPL_H

#include <functional>
#include <vector>

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

} // namespace gmx

#endif
