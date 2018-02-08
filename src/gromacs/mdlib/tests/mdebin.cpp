/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/mdlib/mdebin.h"

#include <cstdio>

#include <gtest/gtest.h>

#include "gromacs/utility/textreader.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! Wraps fclose to discard the return value to use it as a deleter with gmx::unique_cptr.
void fcloseWrapper(FILE *fp)
{
    fclose(fp);
}

class PrintEbin : public ::testing::Test
{
    public:
        t_ebin                         ebin_;
        unique_cptr<t_ebin, done_ebin> ebinGuard_;
        t_mdebin                       mdebin_;

        TestFileManager                fileManager_;
        // TODO This will be more elegant (and run faster) when we
        // refactor the output routines to write to a stream
        // interface, which can already be handled in-memory when
        // running tests.
        std::string logFilename_;
        FILE       *log_;
        unique_cptr<FILE, fcloseWrapper> logFileGuard_;

        TestReferenceData                refData_;
        TestReferenceChecker             checker_;

        PrintEbin() : ebin_ {0}, ebinGuard_(&ebin_), mdebin_ {0},
        fileManager_(), logFilename_(fileManager_.getTemporaryFilePath(".log")),
        log_(std::fopen(logFilename_.c_str(), "w")), logFileGuard_(log_),
        refData_ {}, checker_(refData_.rootChecker())
        {
            mdebin_.ebin = &ebin_;
        }
};

TEST_F(PrintEbin, print_ebin_HandlesEmptyAverages)
{
    ASSERT_NE(log_, nullptr);

    // Test print_ebin
    print_ebin(nullptr, false, false, false, log_,
               0, 0, eprAVER,
               &mdebin_, nullptr, nullptr, nullptr, nullptr);

    // We need to close the file before the contents are available.
    logFileGuard_.reset(nullptr);

    checker_.checkString(TextReader::readFileToString(logFilename_), "log");
}

// TODO Ideally we'd test for proper output of a normal t_mdebin, but
// building one is a lot of work for limited value.

TEST_F(PrintEbin, pr_ebin_HandlesAverages)
{
    ASSERT_NE(log_, nullptr);

    // Set up the energy entries
    const char *firstName[]      = {"first"};
    const char *secondName[]     = {"second"};
    int         first            = get_ebin_space(mdebin_.ebin, 1, firstName, nullptr);
    int         second           = get_ebin_space(mdebin_.ebin, 1, secondName, nullptr);

    // Put some data into the energy entries
    const real  timevalues[2][2] = { { 1.0, 10.0}, {2.0, 20.0} };
    gmx_bool    bSum             = true;
    for (const auto &values : timevalues)
    {
        add_ebin(mdebin_.ebin, first, 1, &values[0], bSum);
        add_ebin(mdebin_.ebin, second, 1, &values[1], bSum);
        ebin_increase_count(mdebin_.ebin, bSum);
    }

    // Test pr_ebin
    pr_ebin(log_, &ebin_, 0, 2, 5, eprAVER, true);

    // We need to close the file before the contents are available.
    logFileGuard_.reset(nullptr);

    checker_.checkString(TextReader::readFileToString(logFilename_), "log");
}

TEST_F(PrintEbin, pr_ebin_HandlesEmptyAverages)
{
    ASSERT_NE(log_, nullptr);

    // Set up the energy entries
    const char *firstName[]      = {"first"};
    const char *secondName[]     = {"second"};
    get_ebin_space(mdebin_.ebin, 1, firstName, nullptr);
    get_ebin_space(mdebin_.ebin, 1, secondName, nullptr);

    // Test pr_ebin
    pr_ebin(log_, &ebin_, 0, 2, 5, eprAVER, true);

    // We need to close the file before the contents are available.
    logFileGuard_.reset(nullptr);

    checker_.checkString(TextReader::readFileToString(logFilename_), "log");
}

}      // namespace
}      // namespace
}      // namespace
