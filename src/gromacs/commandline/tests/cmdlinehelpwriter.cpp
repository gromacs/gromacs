/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
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
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_commandline
 */
#include <gtest/gtest.h>

#include "gromacs/legacyheaders/types/simple.h"

#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionfileoption.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/utility/file.h"

#include "testutils/stringtest.h"
#include "testutils/testfilemanager.h"

namespace
{

class CommandLineHelpWriterTest : public ::gmx::test::StringTestBase
{
    public:
        void checkHelp(gmx::CommandLineHelpWriter *writer);

        gmx::test::TestFileManager tempFiles_;
};

void CommandLineHelpWriterTest::checkHelp(gmx::CommandLineHelpWriter *writer)
{
    std::string filename = tempFiles_.getTemporaryFilePath("helptext.txt");
    gmx::File file(filename, "w");
    gmx::HelpWriterContext context(&file, gmx::eHelpOutputFormat_Console);
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
    const char * const enumValues[] = {"no", "opt1", "opt2", NULL};
    options.addOption(StringOption("enum").description("Enum option")
                        .enumValue(enumValues).defaultEnumIndex(0));

    std::string filename;
    options.addOption(FileNameOption("f")
                        .description("Input file description")
                        .filetype(eftTrajectory).inputFile().required()
                        .defaultValue("traj"));
    options.addOption(FileNameOption("lib")
                        .description("Library file description")
                        .filetype(eftGenericData).inputFile().libraryFile()
                        .defaultValueIfSet("libdata"));
    options.addOption(FileNameOption("io")
                        .store(&filename)
                        .description("Input/Output file description")
                        .filetype(eftGenericData).inputOutputFile()
                        .defaultValueIfSet("inout"));
    options.addOption(FileNameOption("o")
                        .description("Output file description")
                        .filetype(eftPlot).outputFile());

    options.addOption(SelectionFileOption("sf"));
    options.addOption(SelectionOption("sel").description("Selection option"));

    CommandLineHelpWriter writer(options);
    writer.setShowHidden(true);
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
                        .defaultValue("path/to/long/trajectory/name"));
    options.addOption(FileNameOption("f2")
                        .description("File name option with a long value")
                        .filetype(eftTrajectory).inputFile().required()
                        .defaultValue("path/to/long/trajectory"));
    options.addOption(FileNameOption("lib")
                        .description("File name option with a long value and type")
                        .filetype(eftTrajectory).inputFile().libraryFile()
                        .defaultValue("path/to/long/trajectory/name"));
    options.addOption(FileNameOption("longfileopt")
                        .description("File name option with a long name")
                        .filetype(eftGenericData).inputFile()
                        .defaultValue("deffile"));
    options.addOption(FileNameOption("longfileopt2")
                        .description("File name option with multiple long fields")
                        .filetype(eftGenericData).inputFile().libraryFile()
                        .defaultValue("path/to/long/file/name"));

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

/*
 * Tests help printing with selection options with values.
 */
TEST_F(CommandLineHelpWriterTest, HandlesSelectionOptions)
{
    using gmx::SelectionFileOption;
    using gmx::SelectionOption;

    gmx::Options options(NULL, NULL);
    options.addOption(SelectionFileOption("sf"));
    options.addOption(SelectionOption("refsel").required()
                        .description("Reference selection option"));
    options.addOption(SelectionOption("sel").required().valueCount(2)
                        .description("Selection option"));
    gmx::SelectionCollection selections;
    gmx::SelectionOptionManager manager(&selections);
    setManagerForSelectionOptions(&options, &manager);
    options.finish();
    manager.parseRequestedFromString(
            "resname SOL;"
            "surface = within 0.5 of resname SOL;"
            "group \"Protein\" and surface;"
            "group \"Protein\" and not surface;");

    gmx::CommandLineHelpWriter writer(options);
    checkHelp(&writer);
}

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
