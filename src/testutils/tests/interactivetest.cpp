/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Self-tests for interactive test helpers.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/interactivetest.h"

#include <cstddef>

#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/textstream.h"

#include "testutils/refdata.h"

namespace gmx
{
namespace test
{
namespace
{

class InteractiveSession
{
public:
    explicit InteractiveSession(ReferenceDataMode mode) :
        data_(mode), helper_(data_.rootChecker()), nextInputLine_(0)
    {
    }

    void addOutput(const char* output) { events_.emplace_back(WriteOutput, output); }
    void addInputLine(const char* inputLine) { inputLines_.push_back(inputLine); }
    void addReadInput() { events_.emplace_back(ReadInput, ""); }
    void addInput(const char* inputLine)
    {
        addInputLine(inputLine);
        addReadInput();
    }
    void addInputNoNewline(const char* inputLine)
    {
        addInputLine(inputLine);
        helper_.setLastNewline(false);
        events_.emplace_back(ReadInputNoNewline, "");
    }

    void run()
    {
        gmx::TextInputStream&  input  = helper_.inputStream();
        gmx::TextOutputStream& output = helper_.outputStream();
        helper_.setInputLines(inputLines_);
        std::vector<Event>::const_iterator event;
        for (event = events_.begin(); event != events_.end(); ++event)
        {
            if (event->first == WriteOutput)
            {
                output.write(event->second);
            }
            else
            {
                std::string expectedLine;
                const bool  bInputRemaining = (nextInputLine_ < inputLines_.size());
                if (bInputRemaining)
                {
                    expectedLine = inputLines_[nextInputLine_];
                    if (event->first != ReadInputNoNewline)
                    {
                        expectedLine.append("\n");
                    }
                }
                ++nextInputLine_;
                std::string line;
                EXPECT_EQ(bInputRemaining, input.readLine(&line));
                EXPECT_EQ(expectedLine, line);
            }
        }
        helper_.checkSession();
    }

private:
    enum EventType
    {
        ReadInput,
        ReadInputNoNewline,
        WriteOutput
    };
    // The latter is the output string.
    typedef std::pair<EventType, const char*> Event;

    TestReferenceData        data_;
    InteractiveTestHelper    helper_;
    std::vector<const char*> inputLines_;
    size_t                   nextInputLine_;
    std::vector<Event>       events_;
};

TEST(InteractiveTestHelperTest, ChecksSimpleSession)
{
    {
        InteractiveSession session(ReferenceDataMode::UpdateAll);
        session.addOutput("First line\n");
        session.addOutput("> ");
        session.addInput("input");
        session.addOutput("Second line\n");
        session.addOutput("> ");
        session.addReadInput();
        session.addOutput("\n");
        session.addOutput(".\n");
        session.run();
    }
    {
        InteractiveSession session(ReferenceDataMode::Compare);
        session.addOutput("First line\n");
        session.addOutput("> ");
        session.addInput("input");
        session.addOutput("Second line\n");
        session.addOutput("> ");
        session.addReadInput();
        session.addOutput("\n");
        session.addOutput(".\n");
        session.run();
    }
}

TEST(InteractiveTestHelperTest, ChecksSessionWithoutLastNewline)
{
    {
        InteractiveSession session(ReferenceDataMode::UpdateAll);
        session.addOutput("First line\n");
        session.addOutput("> ");
        session.addInput("input");
        session.addOutput("Second line\n");
        session.addOutput("> ");
        session.addInputNoNewline("input2");
        session.addOutput("\n");
        session.addOutput(".\n");
        session.run();
    }
    {
        InteractiveSession session(ReferenceDataMode::Compare);
        session.addOutput("First line\n");
        session.addOutput("> ");
        session.addInput("input");
        session.addOutput("Second line\n");
        session.addOutput("> ");
        session.addInputNoNewline("input2");
        session.addOutput("\n");
        session.addOutput(".\n");
        session.run();
    }
}

TEST(InteractiveTestHelperTest, ChecksSessionWithMissingOutput)
{
    {
        InteractiveSession session(ReferenceDataMode::UpdateAll);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addInput("input2");
        session.addOutput("Second line\n> ");
        session.addReadInput();
        session.addOutput("\n.\n");
        session.run();
    }
    {
        InteractiveSession session(ReferenceDataMode::Compare);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addInput("input2");
        session.addOutput("Second line\n> ");
        session.addReadInput();
        session.addOutput("\n.\n");
        session.run();
    }
}

TEST(InteractiveTestHelperTest, ChecksSessionWithEquivalentOutput)
{
    {
        InteractiveSession session(ReferenceDataMode::UpdateAll);
        session.addOutput("First line\n");
        session.addOutput("> ");
        session.addInput("input");
        session.addOutput("Second line\n> ");
        session.addReadInput();
        session.addOutput("\n");
        session.addOutput(".\n");
        session.run();
    }
    {
        InteractiveSession session(ReferenceDataMode::Compare);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addOutput("Second line\n");
        session.addOutput("> ");
        session.addReadInput();
        session.addOutput("\n.\n");
        session.run();
    }
}

TEST(InteractiveTestHelperTest, DetectsIncorrectOutput)
{
    {
        InteractiveSession session(ReferenceDataMode::UpdateAll);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addOutput("Second line\n> ");
        session.addReadInput();
        session.addOutput("\n.\n");
        session.run();
    }
    {
        InteractiveSession session(ReferenceDataMode::Compare);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addOutput("Incorrect line\n> ");
        session.addReadInput();
        session.addOutput("\n.\n");
        EXPECT_NONFATAL_FAILURE(session.run(), "");
    }
}

TEST(InteractiveTestHelperTest, DetectsMissingOutput)
{
    {
        InteractiveSession session(ReferenceDataMode::UpdateAll);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addOutput("Second line\n> ");
        session.addInput("input2");
        session.addOutput("Third line\n> ");
        session.addReadInput();
        session.addOutput("\n.\n");
        session.run();
    }
    {
        InteractiveSession session(ReferenceDataMode::Compare);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addInput("input2");
        session.addOutput("Third line\n> ");
        session.addReadInput();
        session.addOutput("\n.\n");
        EXPECT_NONFATAL_FAILURE(session.run(), "");
    }
}

TEST(InteractiveTestHelperTest, DetectsMissingFinalOutput)
{
    {
        InteractiveSession session(ReferenceDataMode::UpdateAll);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addOutput("Second line\n> ");
        session.addReadInput();
        session.addOutput("\n.\n");
        session.run();
    }
    {
        InteractiveSession session(ReferenceDataMode::Compare);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addOutput("Second line\n> ");
        session.addReadInput();
        EXPECT_NONFATAL_FAILURE(session.run(), "");
    }
}

TEST(InteractiveTestHelperTest, DetectsExtraOutput)
{
    {
        InteractiveSession session(ReferenceDataMode::UpdateAll);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addInput("input2");
        session.addOutput("More output\n> ");
        session.addReadInput();
        session.addOutput("\n.\n");
        session.run();
    }
    {
        InteractiveSession session(ReferenceDataMode::Compare);
        session.addOutput("First line\n> ");
        session.addInput("input");
        session.addOutput("Extra output\n> ");
        session.addInput("input2");
        session.addOutput("More output\n> ");
        session.addReadInput();
        session.addOutput("\n.\n");
        EXPECT_NONFATAL_FAILURE(session.run(), "");
    }
}

TEST(InteractiveTestHelperTest, DetectsMissingInput)
{
    {
        InteractiveSession session(ReferenceDataMode::UpdateAll);
        session.addInput("input");
        session.addInput("input2");
        session.addReadInput();
        session.run();
    }
    {
        InteractiveSession session(ReferenceDataMode::Compare);
        session.addInputLine("input");
        session.addInputLine("input2");
        session.addReadInput();
        session.addReadInput();
        EXPECT_NONFATAL_FAILURE(session.run(), "");
    }
}

TEST(InteractiveTestHelperTest, DetectsExtraInput)
{
    {
        InteractiveSession session(ReferenceDataMode::UpdateAll);
        session.addInput("input");
        session.addInput("input2");
        session.addReadInput();
        session.run();
    }
    {
        InteractiveSession session(ReferenceDataMode::Compare);
        session.addInputLine("input");
        session.addInputLine("input2");
        session.addReadInput();
        session.addReadInput();
        session.addReadInput();
        session.addReadInput();
        EXPECT_NONFATAL_FAILURE(session.run(), "");
    }
}

} // namespace
} // namespace test
} // namespace gmx
