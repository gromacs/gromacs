/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2019, by the GROMACS development team, led by
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
 * Tests gmx::CommandLineHelpWriter.
 *
 * These tests fail for any change in the output, and it should be reviewed
 * whether the change was intentional.
 * For development, the tests can be run with a '-stdout' command-line option
 * to print out the help to stdout instead of using the XML reference
 * framework.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "gromacs/commandline/cmdlinehelpwriter.h"

#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/stringtest.h"

namespace
{

class CommandLineHelpWriterTest : public ::gmx::test::StringTestBase
{
public:
    CommandLineHelpWriterTest() : bHidden_(false) {}

    void checkHelp(gmx::CommandLineHelpWriter* writer);

    bool bHidden_;
};

void CommandLineHelpWriterTest::checkHelp(gmx::CommandLineHelpWriter* writer)
{
    gmx::StringOutputStream     stream;
    gmx::TextWriter             streamWriter(&stream);
    gmx::CommandLineHelpContext context(&streamWriter, gmx::eHelpOutputFormat_Console, nullptr,
                                        "test");
    context.setShowHidden(bHidden_);
    writer->writeHelp(context);
    stream.close();

    checkText(stream.toString(), "HelpText");
}


/********************************************************************
 * Tests start here
 */

/*
 * Tests help printing for each option type, but doesn't contain much
 * variablity in the options.
 */
TEST_F(CommandLineHelpWriterTest, HandlesOptionTypes)
{
    using namespace gmx;

    Options options;
    options.addOption(BooleanOption("bool").description("Boolean option").defaultValue(true));
    options.addOption(BooleanOption("hidden").description("Hidden option").hidden().defaultValue(true));
    options.addOption(IntegerOption("int").description("Integer option").defaultValue(2));
    ivec intvec = { 1, 2, 3 };
    options.addOption(IntegerOption("ivec").description("Integer vector option").vector().store(intvec));
    options.addOption(DoubleOption("double").description("Double option").defaultValue(2.5));
    dvec dblvec = { 1.1, 2.3, 3.2 };
    options.addOption(DoubleOption("dvec").description("Double vector option").vector().store(dblvec));
    options.addOption(DoubleOption("time").description("Time option (%t)").timeValue().defaultValue(10.0));
    options.addOption(StringOption("string").description("String option").defaultValue("test"));
    const char* const enumValues[] = { "no", "opt1", "opt2" };
    options.addOption(
            StringOption("enum").description("Enum option").enumValue(enumValues).defaultEnumIndex(0));
    options.addOption(
            EnumIntOption("ienum").description("Enum option").enumValue(enumValues).defaultValue(1));

    std::string filename;
    options.addOption(FileNameOption("f")
                              .description("Input file description")
                              .filetype(eftTrajectory)
                              .inputFile()
                              .required()
                              .defaultBasename("traj"));
    options.addOption(FileNameOption("mult")
                              .description("Multiple file description")
                              .filetype(eftTrajectory)
                              .inputFile()
                              .multiValue()
                              .defaultBasename("traj"));
    options.addOption(FileNameOption("lib")
                              .description("Library file description")
                              .filetype(eftGenericData)
                              .inputFile()
                              .libraryFile()
                              .defaultBasename("libdata"));
    options.addOption(FileNameOption("io")
                              .store(&filename)
                              .description("Input/Output file description")
                              .filetype(eftGenericData)
                              .inputOutputFile()
                              .defaultBasename("inout"));
    options.addOption(
            FileNameOption("o").description("Output file description").filetype(eftPlot).outputFile());

    CommandLineHelpWriter writer(options);
    bHidden_ = true;
    checkHelp(&writer);
}

//! Enum value for testing.
enum TestEnum
{
    eFoo,
    eBar
};

/*
 * Tests that default values taken from variables are properly visible in the
 * help output.
 */
TEST_F(CommandLineHelpWriterTest, HandlesDefaultValuesFromVariables)
{
    using namespace gmx;

    Options options;

    bool bValue = true;
    options.addOption(BooleanOption("bool").description("Boolean option").store(&bValue));

    int ivalue = 3;
    options.addOption(IntegerOption("int").description("Integer option").store(&ivalue));

    int iavalue[] = { 2, 3 };
    options.addOption(
            IntegerOption("int2").description("Integer 2-value option").store(iavalue).valueCount(2));

    std::vector<std::string> svalues;
    svalues.emplace_back("foo");
    options.addOption(StringOption("str").description("String option").storeVector(&svalues).multiValue());

    TestEnum          evalue    = eBar;
    const char* const allowed[] = { "foo", "bar" };
    options.addOption(
            EnumOption<TestEnum>("enum").description("Enum option").enumValue(allowed).store(&evalue));

    CommandLineHelpWriter writer(options);
    checkHelp(&writer);
}

/*
 * Tests help printing with file name options with various values that don't
 * fit into the allocated columns.
 */
TEST_F(CommandLineHelpWriterTest, HandlesLongFileOptions)
{
    using gmx::eftGenericData;
    using gmx::eftTrajectory;
    using gmx::FileNameOption;

    gmx::Options options;
    options.addOption(FileNameOption("f")
                              .description("File name option with a long value")
                              .filetype(eftTrajectory)
                              .inputFile()
                              .required()
                              .defaultBasename("path/to/long/trajectory/name"));
    options.addOption(FileNameOption("f2")
                              .description("File name option with a long value")
                              .filetype(eftTrajectory)
                              .inputFile()
                              .required()
                              .defaultBasename("path/to/long/trajectory"));
    options.addOption(FileNameOption("lib")
                              .description("File name option with a long value and type")
                              .filetype(eftTrajectory)
                              .inputFile()
                              .libraryFile()
                              .defaultBasename("path/to/long/trajectory/name"));
    options.addOption(FileNameOption("longfileopt")
                              .description("File name option with a long name")
                              .filetype(eftGenericData)
                              .inputFile()
                              .defaultBasename("deffile"));
    options.addOption(FileNameOption("longfileopt2")
                              .description("File name option with multiple long fields")
                              .filetype(eftGenericData)
                              .inputFile()
                              .libraryFile()
                              .defaultBasename("path/to/long/file/name"));

    gmx::CommandLineHelpWriter writer(options);
    checkHelp(&writer);
}

/*
 * Tests help printing with general options with various values that don't
 * fit into the allocated columns.
 */
TEST_F(CommandLineHelpWriterTest, HandlesLongOptions)
{
    using gmx::BooleanOption;
    using gmx::DoubleOption;
    using gmx::StringOption;

    gmx::Options options;
    options.addOption(
            BooleanOption("longboolean").description("Boolean option with a long name").defaultValue(true));
    dvec dblvec = { 1.135, 2.32, 3.2132 };
    options.addOption(DoubleOption("dvec").description("Double vector option").vector().store(dblvec));
    std::vector<std::string> values;
    values.emplace_back("A very long string value that overflows even the description column");
    values.emplace_back(
            "Another very long string value that overflows even the description column");
    options.addOption(StringOption("string")
                              .description("String option with very long values (may "
                                           "be less relevant with selections having "
                                           "their own option type)")
                              .storeVector(&values));

    gmx::CommandLineHelpWriter writer(options);
    checkHelp(&writer);
}

/* TODO: Add corresponding tests to either the selection module, or as part of
 * trajectoryanalysis tests.
 * Tests help printing with selection options with values.
 */
#if 0
TEST_F(CommandLineHelpWriterTest, HandlesSelectionOptions)
{
    using gmx::SelectionFileOption;
    using gmx::SelectionOption;

    gmx::Options                options;
    gmx::SelectionCollection    selections;
    gmx::SelectionOptionManager manager(&selections);
    options.addManager(&manager);
    options.addOption(SelectionFileOption("sf"));
    options.addOption(SelectionOption("refsel").required()
                          .description("Reference selection option"));
    options.addOption(SelectionOption("sel").required().valueCount(2)
                          .description("Selection option"));
    options.finish();
    manager.parseRequestedFromString(
            "resname SOL;"
            "surface = within 0.5 of resname SOL;"
            "group \"Protein\" and surface;"
            "group \"Protein\" and not surface;");

    gmx::CommandLineHelpWriter writer(options);
    checkHelp(&writer);
}
#endif

/*
 * Tests help output with option groups.
 */
TEST_F(CommandLineHelpWriterTest, HandlesOptionGroups)
{
    using gmx::IntegerOption;

    gmx::Options            options;
    gmx::IOptionsContainer& group1 = options.addGroup();
    gmx::IOptionsContainer& group2 = options.addGroup();
    group2.addOption(IntegerOption("sub2").description("Option in group 2"));
    group1.addOption(IntegerOption("sub11").description("Option in group 1"));
    options.addOption(IntegerOption("main").description("Option in root group"));
    group1.addOption(IntegerOption("sub12").description("Option in group 1"));

    gmx::CommandLineHelpWriter writer(options);
    checkHelp(&writer);
}

/*
 * Tests help output using a help text.
 */
TEST_F(CommandLineHelpWriterTest, HandlesHelpText)
{
    const char* const help[] = { "Help text", "for testing." };
    using gmx::IntegerOption;

    gmx::Options options;
    options.addOption(IntegerOption("int").description("Integer option").defaultValue(2));

    gmx::CommandLineHelpWriter writer(options);
    writer.setHelpText(help);
    checkHelp(&writer);
}

/*
 * Test known issue output.
 */
TEST_F(CommandLineHelpWriterTest, HandlesKnownIssues)
{
    const char* const bugs[] = { "This is a bug.", "And this is another one." };
    using gmx::IntegerOption;

    gmx::Options options;
    options.addOption(IntegerOption("int").description("Integer option").defaultValue(2));

    gmx::CommandLineHelpWriter writer(options);
    writer.setKnownIssues(bugs);
    checkHelp(&writer);
}

} // namespace
