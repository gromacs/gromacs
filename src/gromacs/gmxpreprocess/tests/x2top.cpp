/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests for gmx x2top.
 *
 * The tests run x2top on small structures and compare the topology that is written with
 * reference data. They cover the default OPLS-AA behavior (which has to stay as it is), and
 * the options for MINFF: the compact MINFF format, -n2t, -ktH, -dH, -dM, -aH and -smear.
 * The MINFF tests use a small min.ff directory with the name-to-type table in the test data.
 *
 * \author Michael Holmboe <michael.holmboe@umu.se>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/x2top.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Strings containing regular expressions for lines to skip when matching.
 *
 * These are the header lines of the topology that depend on the machine, the time and the
 * location of the files. */
const std::vector<std::string> c_regexStringsToSkip = { "^;[[:blank:]] *File '.*' was generated.*",
                                                        "^;[[:blank:]]*By user:.*",
                                                        "^;[[:blank:]]*On host:.*",
                                                        "^;[[:blank:]]*At date:.*",
                                                        "^;[[:blank:]]*:-\\).*\\(-:.*",
                                                        "^;[[:blank:]]*Executable:.*",
                                                        "^;[[:blank:]]*Data prefix:.*",
                                                        "^;[[:blank:]]*Working dir:.*",
                                                        "^;[[:blank:]]*Command line:.*",
                                                        "^;[[:blank:]]*[[:alnum:]_-]*-test.*" };

//! Compiled regular expressions for lines to skip when matching.
FilteringExactTextMatch c_textMatcher(c_regexStringsToSkip, false, true);

/*! \brief Makes a directory the working directory for the lifetime of the object.
 *
 * x2top finds force fields in the working directory first, which lets the tests use the
 * small min.ff of the test data. All files that the tests use are given with full paths. */
class ScopedWorkingDirectory
{
public:
    explicit ScopedWorkingDirectory(const std::filesystem::path& directory) :
        previous_(std::filesystem::current_path())
    {
        std::filesystem::current_path(directory);
    }
    ~ScopedWorkingDirectory() { std::filesystem::current_path(previous_); }

    ScopedWorkingDirectory(const ScopedWorkingDirectory&)            = delete;
    ScopedWorkingDirectory& operator=(const ScopedWorkingDirectory&) = delete;

private:
    std::filesystem::path previous_;
};

class X2topTest : public CommandLineTestBase
{
public:
    X2topTest() { setOutputFile("-o", "topol.top", TextFileMatch(c_textMatcher)); }

    /*! \brief Runs x2top on a structure from the test data.
     *
     * \param[in] structure   Name of the structure file in the test data directory.
     * \param[in] options     Further command-line options.
     * \param[in] useMinFF    Whether to use the min.ff in the test data, by running in its directory.
     */
    void runTest(const std::string& structure, const std::vector<std::string>& options, const bool useMinFF)
    {
        setInputFile("-f", ("x2top_data/" + structure).c_str());

        CommandLine& cmdline = commandLine();
        // Provide the name of the module to call
        std::string module[] = { "x2top" };
        cmdline.merge(CommandLine(module));
        cmdline.merge(CommandLine(options));

        std::unique_ptr<ScopedWorkingDirectory> workingDirectory;
        if (useMinFF)
        {
            workingDirectory = std::make_unique<ScopedWorkingDirectory>(
                    TestFileManager::getInputDataDirectory() / "x2top_data");
        }
        ASSERT_EQ(0, gmx_x2top(cmdline.argc(), cmdline.argv()));
        checkOutputFiles();
    }

    //! Returns the full path of a file in the x2top test data, to give to an option.
    static std::string dataFile(const std::string& name)
    {
        return (TestFileManager::getInputDataDirectory() / "x2top_data" / name).string();
    }
};

// Behavior that has to stay as it is for force fields other than MINFF

TEST_F(X2topTest, OplsaaEthane)
{
    runTest("ethane.pdb", { "-ff", "oplsaa", "-name", "ETH" }, false);
}

TEST_F(X2topTest, OplsaaEthaneWrappedAcrossBoundary)
{
    // The molecule is split over the periodic boundary, so this checks the bonds and angles
    // that are found across it.
    runTest("ethane_wrapped.pdb", { "-ff", "oplsaa", "-name", "ETH" }, false);
}

TEST_F(X2topTest, OplsaaEthaneWithHydrogenAngleConstant)
{
    runTest("ethane.pdb", { "-ff", "oplsaa", "-name", "ETH", "-kt", "300", "-ktH", "123" }, false);
}

// Compact MINFF format

TEST_F(X2topTest, MinffDefaultFormat)
{
    runTest("kaolinite.pdb", { "-ff", "min" }, true);
}

TEST_F(X2topTest, MinffAliasWithBondConstant)
{
    runTest("kaolinite.pdb", { "-ff", "minff", "-kb", "441050" }, true);
}

TEST_F(X2topTest, MinffHydrogenBondAndAngleOverrides)
{
    runTest("kaolinite.pdb", { "-ff", "min", "-dH", "0.0957", "-aH", "109.47" }, true);
}

TEST_F(X2topTest, MinffBondLengthOverrideAndAngleConstants)
{
    runTest("kaolinite.pdb", { "-ff", "min", "-dM", "0.15", "-kt", "450", "-ktH", "90" }, true);
}

TEST_F(X2topTest, MinffHydrogenFoundForAtomNamedHo)
{
    // The hydroxyl hydrogens are named Ho (not holmium) and have Ho as element in the PDB
    // file, so they have to be found from the mass in the name-to-type table.
    runTest("kaolinite_Ho.pdb", { "-ff", "min", "-aH", "109.5" }, true);
}

TEST_F(X2topTest, MinffHydrogenFoundForAtomNamedHoInGroFile)
{
    // A gro file has no element information
    runTest("kaolinite_Ho.gro", { "-ff", "min", "-aH", "109.5" }, true);
}

// Name-to-type table from a file, and charge smearing

TEST_F(X2topTest, MinffTableFromFile)
{
    // The oxygen charges of this table are not right for the structure
    runTest("kaolinite.pdb", { "-ff", "min", "-n2t", dataFile("kaolinite_dummyO.n2t") }, true);
}

TEST_F(X2topTest, MinffSmearedChargesFromTableWithWrongOxygenCharges)
{
    // The total charge is not an integer without -smear, see the previous test
    runTest("kaolinite.pdb", { "-ff", "min", "-n2t", dataFile("kaolinite_dummyO.n2t"), "-smear" }, true);
}

TEST_F(X2topTest, MinffSmearedChargesFromTableThatIsRight)
{
    // The charges of the table are the smeared ones, so they stay as they are
    runTest("kaolinite.pdb", { "-ff", "min", "-smear" }, true);
}

} // namespace
} // namespace test
} // namespace gmx
