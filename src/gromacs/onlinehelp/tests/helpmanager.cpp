/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Tests for help topic management and help topic formatting.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#include "gmxpre.h"

#include "gromacs/onlinehelp/helpmanager.h"

#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/onlinehelp/helptopic.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"

#include "gromacs/onlinehelp/tests/mock_helptopic.h"
#include "testutils/stringtest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace
{

using gmx::test::MockHelpTopic;

class HelpTestBase : public gmx::test::StringTestBase
{
    public:
        HelpTestBase();

        gmx::test::TestFileManager tempFiles_;
        MockHelpTopic              rootTopic_;
        std::string                filename_;
        gmx::File                  helpFile_;
        gmx::HelpWriterContext     context_;
        gmx::HelpManager           manager_;
};

HelpTestBase::HelpTestBase()
    : rootTopic_("", NULL, "Root topic text"),
      filename_(tempFiles_.getTemporaryFilePath("helptext.txt")),
      helpFile_(filename_, "w"),
      context_(&helpFile_, gmx::eHelpOutputFormat_Console),
      manager_(rootTopic_, context_)
{
}

/********************************************************************
 * Tests for HelpManager
 */

//! Test fixture for gmx::HelpManager.
typedef HelpTestBase HelpManagerTest;

TEST_F(HelpManagerTest, HandlesRootTopic)
{
    using ::testing::_;
    EXPECT_CALL(rootTopic_, writeHelp(_));
    manager_.writeCurrentTopic();
}

TEST_F(HelpManagerTest, HandlesSubTopics)
{
    MockHelpTopic &first =
        rootTopic_.addSubTopic("first", "First topic", NULL);
    MockHelpTopic &firstSub =
        first.addSubTopic("firstsub", "First subtopic", NULL);
    rootTopic_.addSubTopic("second", "Second topic", NULL);

    using ::testing::_;
    EXPECT_CALL(firstSub, writeHelp(_));
    ASSERT_NO_THROW_GMX(manager_.enterTopic("first"));
    ASSERT_NO_THROW_GMX(manager_.enterTopic("firstsub"));
    manager_.writeCurrentTopic();
}

TEST_F(HelpManagerTest, HandlesInvalidTopics)
{
    MockHelpTopic &first =
        rootTopic_.addSubTopic("first", "First topic", NULL);
    first.addSubTopic("firstsub", "First subtopic", NULL);
    rootTopic_.addSubTopic("second", "Second topic", NULL);

    ASSERT_THROW_GMX(manager_.enterTopic("unknown"), gmx::InvalidInputError);
    ASSERT_NO_THROW_GMX(manager_.enterTopic("first"));
    ASSERT_THROW_GMX(manager_.enterTopic("unknown"), gmx::InvalidInputError);
    ASSERT_THROW_GMX(manager_.enterTopic("second"), gmx::InvalidInputError);
    ASSERT_NO_THROW_GMX(manager_.enterTopic("firstsub"));
}

/********************************************************************
 * Tests for help topic formatting
 */

struct TestHelpText
{
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

const char        TestHelpText::name[]  = "testtopic";
const char        TestHelpText::title[] = "Topic title";
const char *const TestHelpText::text[]  = {
    "Test topic text.[PAR]",
    "Another paragraph of text."
};

class HelpTopicFormattingTest : public HelpTestBase
{
    public:
        void checkHelpFormatting();
};

void HelpTopicFormattingTest::checkHelpFormatting()
{
    ASSERT_NO_THROW_GMX(manager_.enterTopic("testtopic"));
    ASSERT_NO_THROW_GMX(manager_.writeCurrentTopic());
    helpFile_.close();

    checkFileContents(filename_, "HelpText");
}

TEST_F(HelpTopicFormattingTest, FormatsSimpleTopic)
{
    rootTopic_.addSubTopic(gmx::HelpTopicPointer(
                                   new gmx::SimpleHelpTopic<TestHelpText>));
    checkHelpFormatting();
}

TEST_F(HelpTopicFormattingTest, FormatsCompositeTopicWithSubTopics)
{
    gmx::CompositeHelpTopicPointer topic(new gmx::CompositeHelpTopic<TestHelpText>);
    MockHelpTopic::addSubTopic(topic.get(), "subtopic", "First subtopic", "Text");
    MockHelpTopic::addSubTopic(topic.get(), "other", "Second subtopic", "Text");
    rootTopic_.addSubTopic(move(topic));
    checkHelpFormatting();
}

} // namespace
