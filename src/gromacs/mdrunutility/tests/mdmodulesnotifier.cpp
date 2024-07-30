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
/*! \internal \file
 * \brief
 * Tests MDModulesNotifier
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include "gromacs/mdrunutility/mdmodulesnotifier.h"

#include <functional>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace gmx
{
namespace test
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

    bool notifiedEventA() const { return notifiedEventA_; }

private:
    bool notifiedEventA_ = false;
};

class EventBCallee final
{
public:
    void callback(EventB* /* bPointer */) { notifiedEventB_ = true; }

    bool notifiedEventB() const { return notifiedEventB_; }

private:
    bool notifiedEventB_ = false;
};

class EventAandBCallee final
{
public:
    void notify(EventB* /* bPointer */) { notifiedEventB_ = true; }

    void callback(EventA /* a */) { notifiedEventA_ = true; }

    bool notifiedEventB() const { return notifiedEventB_; }
    bool notifiedEventA() const { return notifiedEventA_; }

private:
    bool notifiedEventB_ = false;
    bool notifiedEventA_ = false;
};

TEST(MDModulesNotifierTest, AddConsumer)
{
    BuildMDModulesNotifier<EventA>::type notifier;
    EventACallee                         eventACallee;

    EXPECT_FALSE(eventACallee.notifiedEventA());

    notifier.subscribe([&eventACallee](EventA eventA) { eventACallee.callback(eventA); });
    notifier.notify(EventA{});

    EXPECT_TRUE(eventACallee.notifiedEventA());
}

TEST(MDModulesNotifierTest, AddConsumerWithPointerParameter)
{
    BuildMDModulesNotifier<EventB*>::type notifier;
    EventBCallee                          eventBCallee;

    EXPECT_FALSE(eventBCallee.notifiedEventB());

    notifier.subscribe([&eventBCallee](EventB* eventB) { eventBCallee.callback(eventB); });
    EventB* eventBPointer = nullptr;
    notifier.notify(eventBPointer);

    EXPECT_TRUE(eventBCallee.notifiedEventB());
}

TEST(MDModulesNotifierTest, AddTwoDifferentConsumers)
{
    BuildMDModulesNotifier<EventA, EventB*>::type notifier;
    EventBCallee                                  eventBCallee;
    EventACallee                                  eventACallee;

    EXPECT_FALSE(eventACallee.notifiedEventA());
    EXPECT_FALSE(eventBCallee.notifiedEventB());

    notifier.subscribe([&eventBCallee](EventB* eventB) { eventBCallee.callback(eventB); });
    notifier.subscribe([&eventACallee](EventA eventA) { eventACallee.callback(eventA); });

    EventB* eventBPointer = nullptr;
    notifier.notify(eventBPointer);

    EXPECT_FALSE(eventACallee.notifiedEventA());
    EXPECT_TRUE(eventBCallee.notifiedEventB());

    notifier.notify(EventA{});

    EXPECT_TRUE(eventACallee.notifiedEventA());
    EXPECT_TRUE(eventBCallee.notifiedEventB());
}

TEST(MDModulesNotifierTest, AddConsumerOfTwoResources)
{
    BuildMDModulesNotifier<EventA, EventB*>::type notifier;

    EventAandBCallee callee;

    EXPECT_FALSE(callee.notifiedEventB());
    EXPECT_FALSE(callee.notifiedEventA());

    // requires a template parameter here, because call is ambiguous otherwise
    notifier.subscribe([&callee](EventA msg) { callee.callback(msg); });
    notifier.subscribe([&callee](EventB* msg) { callee.notify(msg); });

    EventB* eventBp = nullptr;

    notifier.notify(eventBp);

    EXPECT_FALSE(callee.notifiedEventA());
    EXPECT_TRUE(callee.notifiedEventB());

    notifier.notify(EventA{});

    EXPECT_TRUE(callee.notifiedEventA());
    EXPECT_TRUE(callee.notifiedEventB());
}


} // namespace
} // namespace test
} // namespace gmx
