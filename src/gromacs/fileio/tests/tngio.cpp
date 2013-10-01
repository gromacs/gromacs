/*! \internal \file
 * \brief
 * Tests for file I/O routines
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_fileio
 */
#include <gtest/gtest.h>
#include <string>

#include "../tngio.h"
#include "../tngio_for_tools.h"
#include "testutils/testfilemanager.h"
#include "gromacs/utility/path.h"

namespace
{

class TngTest : public ::testing::Test
{
    public:
        TngTest()
        {
        }
        gmx::test::TestFileManager      fileManager_;
};

TEST_F(TngTest, CanOpenTngFile)
{
    tng_trajectory_t tng;
    // TODO Mock or adapt gmx_fatal so that the error path is not
    // fatal. Or work out how to use assertions properly in the
    // general case. Then, cure cancer.
    tng_open(fileManager_.getInputFilePath("spc2-traj.tng").c_str(),
             'r',
             &tng);
    tng_close(&tng);
}

TEST_F(TngTest, CloseBeforeOpenIsNotFatal)
{
    tng_trajectory_t tng = NULL;
    tng_close(&tng);
}

} // namespace
