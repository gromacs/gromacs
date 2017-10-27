/*
 * FDAGetStressTest.cpp
 *
 *  Created on: Apr 2, 2015
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
        std::string const& result,
        std::string const& reference
    )
     : testDirectory(testDirectory),
       cmdline(cmdline),
       result(result),
       reference(reference)
    {}

    std::string testDirectory;
    std::vector<std::string> cmdline;
    std::string result;
    std::string reference;
};

//! Test fixture for FDA
class FDAGetStress : public ::testing::WithParamInterface<TestDataStructure>,
                     public CommandLineTestBase
{};

//! Test body for FDA
TEST_P(FDAGetStress, Basic)
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
    caller.append("gmx_fda fda_get_stress");
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

    ASSERT_FALSE(gmx_fda_get_stress(caller.argc(), caller.argv()));

    const double error_factor = 1.0e4;
    const bool weight_by_magnitude = false;
    const bool ignore_sign = true;

    LogicallyEqualComparer<weight_by_magnitude,ignore_sign> comparer(error_factor);

    // compare atom pairs
    EXPECT_TRUE((equal(TextSplitter(GetParam().reference), TextSplitter(GetParam().result), comparer)));

    gmx_chdir(cwd.c_str());
}

INSTANTIATE_TEST_CASE_P(AllFDAGetStress, FDAGetStress, ::testing::Values(
    TestDataStructure(
        "glycine_trimer",
        {"-ipf", "fda.pfr"},
        "result.psr",
        "punctual_stress.psr"
    )
));

} // namespace
} // namespace test
} // namespace gmx
