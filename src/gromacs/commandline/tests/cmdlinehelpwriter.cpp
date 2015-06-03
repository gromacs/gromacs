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
#include "gromacs/utility/file.h"

#include "testutils/stringtest.h"
#include "testutils/testfilemanager.h"

namespace
{

class CommandLineHelpWriterTest : public ::gmx::test::StringTestBase
{
    public:
        CommandLineHelpWriterTest() : bHidden_(false) {}

        void checkHelp(gmx::CommandLineHelpWriter *writer);

        gmx::test::TestFileManager tempFiles_;
        bool                       bHidden_;
};

void CommandLineHelpWriterTest::checkHelp(gmx::CommandLineHelpWriter *writer)
{
    std::string                 filename = tempFiles_.getTemporaryFilePath("helptext.txt");
    gmx::File                   file(filename, "w");
    gmx::CommandLineHelpContext context(&file, gmx::eHelpOutputFormat_Console,
                                        NULL, "test");
    context.setShowHidden(bHidden_);
    writer->writeHelp(context);
    file.close();

    checkFileContents(filename, "HelpText");
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

    Options options("test", "Short Description");
    options.addOption(BooleanOption("bool").description("Boolean option")
                          .defaultValue(true));
    options.addOption(BooleanOption("hidden").description("Hidden option")
                          .hidden().defaultValue(true));
    options.addOption(IntegerOption("int").description("Integer option")
                          .defaultValue(2));
    ivec intvec = {1, 2, 3};
    options.addOption(IntegerOption("ivec").description("Integer vector option")
                          .vector().store(intvec));
    options.addOption(DoubleOption("double").description("Double option")
                          .defaultValue(2.5));
    dvec dblvec = {1.1, 2.3, 3.2};
    options.addOption(DoubleOption("dvec").description("Double vector option")
                          .vector().store(dblvec));
    options.addOption(DoubleOption("time").description("Time option (%t)")
                          .timeValue().defaultValue(10.0));
    options.addOption(StringOption("string").description("String option")
                          .defaultValue("test"));
    const char * const enumValues[] = { "no", "opt1", "opt2" };
    options.addOption(StringOption("enum").description("Enum option")
                          .enumValue(enumValues).defaultEnumIndex(0));

    std::string filename;
    options.addOption(FileNameOption("f")
                          .description("Input file description")
                          .filetype(eftTrajectory).inputFile().required()
                          .defaultBasename("traj"));
    options.addOption(FileNameOption("mult")
                          .description("Multiple file description")
                          .filetype(eftTrajectory).inputFile().multiValue()
                          .defaultBasename("traj"));
    options.addOption(FileNameOption("lib")
                          .description("Library file description")
                          .filetype(eftGenericData).inputFile().libraryFile()
                          .defaultBasename("libdata"));
    options.addOption(FileNameOption("io")
                          .store(&filename)
                          .description("Input/Output file description")
                          .filetype(eftGenericData).inputOutputFile()
                          .defaultBasename("inout"));
    options.addOption(FileNameOption("o")
                          .description("Output file description")
                          .filetype(eftPlot).outputFile());

    CommandLineHelpWriter writer(options);
    bHidden_ = true;
    checkHelp(&writer);
}

/*
 * Tests help printing with file name options with various values that don't
 * fit into the allocated columns.
 */
TEST_F(CommandLineHelpWriterTest, HandlesLongFileOptions)
{
    using gmx::FileNameOption;
    using gmx::eftGenericData;
    using gmx::eftTrajectory;

    gmx::Options options(NULL, NULL);
    options.addOption(FileNameOption("f")
                          .description("File name option with a long value")
                          .filetype(eftTrajectory).inputFile().required()
                          .defaultBasename("path/to/long/trajectory/name"));
    options.addOption(FileNameOption("f2")
                          .description("File name option with a long value")
                          .filetype(eftTrajectory).inputFile().required()
                          .defaultBasename("path/to/long/trajectory"));
    options.addOption(FileNameOption("lib")
                          .description("File name option with a long value and type")
                          .filetype(eftTrajectory).inputFile().libraryFile()
                          .defaultBasename("path/to/long/trajectory/name"));
    options.addOption(FileNameOption("longfileopt")
                          .description("File name option with a long name")
                          .filetype(eftGenericData).inputFile()
                          .defaultBasename("deffile"));
    options.addOption(FileNameOption("longfileopt2")
                          .description("File name option with multiple long fields")
                          .filetype(eftGenericData).inputFile().libraryFile()
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

    gmx::Options options(NULL, NULL);
    options.addOption(BooleanOption("longboolean")
                          .description("Boolean option with a long name")
                          .defaultValue(true));
    dvec dblvec = {1.135, 2.32, 3.2132};
    options.addOption(DoubleOption("dvec").description("Double vector option")
                          .vector().store(dblvec));
    std::vector<std::string> values;
    values.push_back("A very long string value that overflows even the description column");
    values.push_back("Another very long string value that overflows even the description column");
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

    gmx::Options                options(NULL, NULL);
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
 * Tests help printing for multiple sections.
 */
TEST_F(CommandLineHelpWriterTest, HandlesMultipleSections)
{
    using namespace gmx;

    Options options("test", "Main Title");
    Options subSect1("subsect1", "Subsection 1 Title");
    Options subSect2("subsect2", "Subsection 2 Title");
    Options subSect3("subsect3", NULL);
    options.addSubSection(&subSect1);
    options.addSubSection(&subSect2);
    options.addSubSection(&subSect3);
    options.setDescription("Description for main section.");
    subSect1.setDescription("Description for subsection 1.");
    subSect2.setDescription("Description for subsection 2.");
    subSect3.setDescription("Description for subsection 3.");
    options.addOption(IntegerOption("main")
                          .description("Option in main section"));
    subSect1.addOption(IntegerOption("sub1")
                           .description("Option in subsection 1"));
    subSect2.addOption(IntegerOption("sub2")
                           .description("Option in subsection 2"));

    CommandLineHelpWriter writer(options);
    writer.setShowDescriptions(true);
    checkHelp(&writer);
}

} // namespace
