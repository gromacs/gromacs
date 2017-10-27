/*
 * FDAViewStressTest.cpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <iostream>
#include <string>
#include <vector>
#include <gtest/gtest.h>
#include "gromacs/gmxana/fda/Helpers.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/path.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "programs/mdrun/mdrun_main.h"
#include "testutils/cmdlinetest.h"
#include "testutils/stdiohelper.h"
#include "testutils/testfilemanager.h"
#include "testutils/TextSplitter.h"
#include "testutils/LogicallyErrorComparer.h"

#define STR(x) #x
#define STRING(x) STR(x)

using namespace fda_analysis;

namespace gmx
{
namespace test
{
namespace
{

struct TestDataStructure
{
    TestDataStructure(
        std::string const& testDirectory,
        std::vector<std::string> const& cmdline,
        std::string const& groupname,
        std::string const& result,
        std::string const& reference
    )
     : testDirectory(testDirectory),
       cmdline(cmdline),
       groupname(groupname),
       result(result),
       reference(reference)
    {}

    std::string testDirectory;
    std::vector<std::string> cmdline;
    std::string groupname;
    std::string result;
    std::string reference;
};

//! Test fixture for FDA
class FDAViewStress : public ::testing::WithParamInterface<TestDataStructure>,
                      public CommandLineTestBase
{};

//! Test body for FDA
TEST_P(FDAViewStress, Basic)
{
    std::string cwd = gmx::Path::getWorkingDirectory();
    std::string dataPath = std::string(STRING(REGRESSIONTEST_PATH)) + "/fda-analysis";
    std::string testPath = fileManager().getTemporaryFilePath("/" + GetParam().testDirectory);

    std::string cmd = "mkdir -p " + testPath;
    ASSERT_FALSE(system(cmd.c_str()));

    cmd = "cp -r " + dataPath + "/" + GetParam().testDirectory + "/* " + testPath;
    ASSERT_FALSE(system(cmd.c_str()));

    gmx_chdir(testPath.c_str());

    ::gmx::test::CommandLine caller;
    caller.append("gmx_fda fda_view_stress");
    for (std::vector<std::string>::const_iterator iterCur(GetParam().cmdline.begin()), iterNext(GetParam().cmdline.begin() + 1),
        iterEnd(GetParam().cmdline.end()); iterCur != iterEnd; ++iterCur, ++iterNext)
    {
        if (iterNext == iterEnd || iterNext->substr(0,1) == "-") caller.append(*iterCur);
        else {
            caller.addOption(iterCur->c_str(), iterNext->c_str());
            ++iterCur, ++iterNext;
        }
    }
    caller.addOption("-o", GetParam().result);

    std::cout << caller.toString() << std::endl;

    if (!GetParam().groupname.empty()) {
        StdioTestHelper stdioHelper(&fileManager());
        stdioHelper.redirectStringToStdin((GetParam().groupname + "\n").c_str());
    }

    ASSERT_FALSE(gmx_fda_view_stress(caller.argc(), caller.argv()));

    const double error_factor = 1.0e4;
    const bool weight_by_magnitude = false;
    const bool ignore_sign = true;

    LogicallyEqualComparer<weight_by_magnitude,ignore_sign> comparer(error_factor);

    // compare atom pairs
    EXPECT_TRUE((equal(TextSplitter(GetParam().reference), TextSplitter(GetParam().result), comparer)));

    gmx_chdir(cwd.c_str());
}

INSTANTIATE_TEST_CASE_P(AllFDAViewStress, FDAViewStress, ::testing::Values(
    TestDataStructure(
        "glycine_trimer",
        {"-f", "punctual_stress.psr", "-frame", "all"},
        "",
        "result.xpm",
        "FDAViewStress.ref0.xpm"
    ),
    TestDataStructure(
        "glycine_trimer",
        {"-f", "virial_stress_von_mises.vma", "-frame", "all"},
        "",
        "result.xpm",
        "FDAViewStress.ref1.xpm"
    ),
    TestDataStructure(
        "glycine_trimer",
        {"-f", "punctual_stress.psa", "-frame", "all"},
        "",
        "result.xpm",
        "FDAViewStress.ref2.xpm"
    ),
    TestDataStructure(
        "glycine_trimer",
        {"-f", "punctual_stress.psa", "-frame", "skip 3"},
        "",
        "result.xpm",
        "FDAViewStress.ref3.xpm"
    ),
    TestDataStructure(
        "glycine_trimer",
        {"-f", "punctual_stress.psa", "-frame", "average 3"},
        "",
        "result.xpm",
        "FDAViewStress.ref4.xpm"
    ),
    TestDataStructure(
        "glycine_trimer",
        {"-f", "punctual_stress.psa", "-s", "glycine_trimer.pdb", "-frame", "0"},
        "",
        "result.pdb",
        "FDAViewStress.ref5.pdb"
    ),
    TestDataStructure(
        "glycine_trimer",
        {"-f", "punctual_stress.psa", "-frame", "3"},
        "",
        "result.xpm",
        "FDAViewStress.ref6.xpm"
    ),
    TestDataStructure(
        "glycine_trimer",
        {"-f", "punctual_stress.psa", "-s", "glycine_trimer.pdb", "-traj", "traj.trr", "-frame", "skip 3"},
        "",
        "result.pdb",
        "FDAViewStress.ref7.pdb"
    )
));

} // namespace
} // namespace test
} // namespace gmx
