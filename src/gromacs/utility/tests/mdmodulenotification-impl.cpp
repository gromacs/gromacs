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
/*! \internal \file
 * \brief
 * Tests MdModuleNotification
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/mdmodulenotification-impl.h"

#include <gmock/gmock.h>

namespace gmx
{

namespace
{

struct EventA
{
};
struct EventB
{
};

class EventACallee final
{
public:
    void callback(EventA /*a*/) { notifiedEventA_ = true; }

    bool notifiedEventA() { return notifiedEventA_; }

private:
    bool notifiedEventA_ = false;
};

class EventBCallee final
{
public:
    void callback(EventB* /* bPointer */) { notifiedEventB_ = true; }

    bool notifiedEventB() { return notifiedEventB_; }

private:
    bool notifiedEventB_ = false;
};

class EventAandBCallee final
{
public:
    void notify(EventB* /* bPointer */) { notifiedEventB_ = true; }

    void callback(EventA /* a */) { notifiedEventA_ = true; }

    bool notifiedEventB() { return notifiedEventB_; }
    bool notifiedEventA() { return notifiedEventA_; }

private:
    bool notifiedEventB_ = false;
    bool notifiedEventA_ = false;
};

TEST(MDModuleNotificationTest, addConsumer)
{
    registerMdModuleNotification<EventA>::type notifications;
    EventACallee                               eventACallee;

    EXPECT_FALSE(eventACallee.notifiedEventA());

    notifications.subscribe([&eventACallee](EventA eventA) { eventACallee.callback(eventA); });
    notifications.notify(EventA{});

    EXPECT_TRUE(eventACallee.notifiedEventA());
}

TEST(MDModuleNotificationTest, addConsumerWithPointerParameter)
{
    registerMdModuleNotification<EventB*>::type notifications;
    EventBCallee                                eventBCallee;

    EXPECT_FALSE(eventBCallee.notifiedEventB());

    notifications.subscribe([&eventBCallee](EventB* eventB) { eventBCallee.callback(eventB); });
    EventB* eventBPointer = nullptr;
    notifications.notify(eventBPointer);

    EXPECT_TRUE(eventBCallee.notifiedEventB());
}

TEST(MDModuleNotificationTest, addTwoDifferentConsumers)
{
    registerMdModuleNotification<EventA, EventB*>::type notifications;
    EventBCallee                                        eventBCallee;
    EventACallee                                        eventACallee;

    EXPECT_FALSE(eventACallee.notifiedEventA());
    EXPECT_FALSE(eventBCallee.notifiedEventB());

    notifications.subscribe([&eventBCallee](EventB* eventB) { eventBCallee.callback(eventB); });
    notifications.subscribe([&eventACallee](EventA eventA) { eventACallee.callback(eventA); });

    EventB* eventBPointer = nullptr;
    notifications.notify(eventBPointer);

    EXPECT_FALSE(eventACallee.notifiedEventA());
    EXPECT_TRUE(eventBCallee.notifiedEventB());

    notifications.notify(EventA{});

    EXPECT_TRUE(eventACallee.notifiedEventA());
    EXPECT_TRUE(eventBCallee.notifiedEventB());
}

TEST(MDModuleNotificationTest, consumerOfTwoResources)
{
    registerMdModuleNotification<EventA, EventB*>::type notifications;

    EventAandBCallee callee;

    EXPECT_FALSE(callee.notifiedEventB());
    EXPECT_FALSE(callee.notifiedEventA());

    // requires a template parameter here, because call is ambiguous otherwise
    notifications.subscribe([&callee](EventA msg) { callee.callback(msg); });
    notifications.subscribe([&callee](EventB* msg) { callee.notify(msg); });

    EventB* eventBp = nullptr;

    notifications.notify(eventBp);

    EXPECT_FALSE(callee.notifiedEventA());
    EXPECT_TRUE(callee.notifiedEventB());

    notifications.notify(EventA{});

    EXPECT_TRUE(callee.notifiedEventA());
    EXPECT_TRUE(callee.notifiedEventB());
}


} // namespace

} // namespace gmx
