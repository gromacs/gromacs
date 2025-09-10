/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief Tests for MessageStringCollector.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/message_string_collector.h"

#include <string>
#include <utility>

#include <gtest/gtest.h>

namespace gmx
{
namespace test
{
namespace
{

TEST(MessageStringCollectorTest, CanAddAndClearMessagesNoContext)
{
    MessageStringCollector messages;

    EXPECT_TRUE(messages.isEmpty());

    messages.append("Message1");

    EXPECT_FALSE(messages.isEmpty());

    messages.append("Message2");

    EXPECT_FALSE(messages.isEmpty());

    messages.clear();

    EXPECT_TRUE(messages.isEmpty());
}

TEST(MessageStringCollectorTest, CanAddAndClearMessagesWithContext)
{
    MessageStringCollector messages;

    messages.startContext("Context1");

    EXPECT_TRUE(messages.isEmpty());

    messages.append("Message1");

    EXPECT_FALSE(messages.isEmpty());

    messages.append("Message1");

    EXPECT_FALSE(messages.isEmpty());

    messages.finishContext();

    messages.clear();

    EXPECT_TRUE(messages.isEmpty());
}

TEST(MessageStringCollectorTest, CanAddStringMessages)
{
    std::string context1 = "Context1";
    std::string message1 = "Message1";

    MessageStringCollector messagesChar;
    MessageStringCollector messagesString;

    messagesChar.startContext(context1.c_str());
    messagesChar.append(message1.c_str());
    messagesChar.finishContext();

    messagesString.startContext(context1);
    messagesString.append(message1);
    messagesString.finishContext();

    EXPECT_EQ(messagesChar.toString(), messagesString.toString());
}

TEST(MessageStringCollectorTest, CanAddCharMessagesConditionally)
{
    std::string context1     = "Context1";
    std::string message1     = "Message1";
    std::string message2     = "Message2";
    bool        conditional1 = true;
    bool        conditional2 = false;

    MessageStringCollector messagesDirect;
    MessageStringCollector messagesConditional;

    messagesDirect.startContext(context1);
    if (conditional1)
    {
        messagesDirect.append(message1.c_str());
    }

    if (conditional2)
    {
        messagesDirect.append(message2.c_str());
    }

    messagesDirect.finishContext();

    messagesConditional.startContext(context1);
    messagesConditional.appendIf(conditional1, message1.c_str());
    messagesConditional.appendIf(conditional2, message2.c_str());
    messagesConditional.finishContext();

    EXPECT_EQ(messagesDirect.toString(), messagesConditional.toString());
}

TEST(MessageStringCollectorTest, CanAddStringMessagesConditionally)
{
    std::string context1     = "Context1";
    std::string message1     = "Message1";
    std::string message2     = "Message2";
    bool        conditional1 = true;
    bool        conditional2 = false;

    MessageStringCollector messagesDirect;
    MessageStringCollector messagesConditional;

    messagesDirect.startContext(context1);
    if (conditional1)
    {
        messagesDirect.append(message1);
    }

    if (conditional2)
    {
        messagesDirect.append(message2);
    }

    messagesDirect.finishContext();

    messagesConditional.startContext(context1);
    messagesConditional.appendIf(conditional1, message1);
    messagesConditional.appendIf(conditional2, message2);
    messagesConditional.finishContext();

    EXPECT_EQ(messagesDirect.toString(), messagesConditional.toString());
}

TEST(MessageStringCollectorTest, CanMoveConstruct)
{
    MessageStringCollector first;
    EXPECT_TRUE(first.isEmpty());
    std::string message = "Message1";
    first.append(message);
    EXPECT_FALSE(first.isEmpty());
    MessageStringCollector second(std::move(first));
    // Now the only valid thing to do with first is to call the
    // destructor.
    EXPECT_FALSE(second.isEmpty());
    EXPECT_EQ(second.toString(), message + "\n");
}

TEST(MessageStringCollectorTest, CanMoveAssign)
{
    MessageStringCollector first, second;
    EXPECT_TRUE(first.isEmpty());
    EXPECT_TRUE(second.isEmpty());
    std::string message = "Message1";
    first.append(message);
    EXPECT_FALSE(first.isEmpty());
    second = std::move(first);
    // Now the only valid thing to do with first is to call the
    // destructor.
    EXPECT_FALSE(second.isEmpty());
    EXPECT_EQ(second.toString(), message + "\n");
}

} // namespace
} // namespace test
} // namespace gmx
