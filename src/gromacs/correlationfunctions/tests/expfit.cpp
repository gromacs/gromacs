/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * Implements test of exponential fitting routines
 *
 * \author Anders G&auml;rden&auml;s <anders.gardenas@gmail.com>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "gromacs/correlationfunctions/expfit.h"

#include "config.h"

#include <cmath>

#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

struct gmx_output_env_t;

namespace gmx
{
namespace test
{
namespace
{
#if HAVE_LMFIT
class ExpfitData
{
public:
    int               nrLines_;
    std::vector<real> x_, y_;
    real              startTime_, endTime_, dt_;
};

class ExpfitTest : public ::testing::Test
{

protected:
    static std::vector<ExpfitData> data_;
    test::TestReferenceData        refData_;
    test::TestReferenceChecker     checker_;
    ExpfitTest() : checker_(refData_.rootChecker()) {}

    // Static initiation, only run once every test.
    static void SetUpTestSuite()
    {
        double**                           tempValues = nullptr;
        std::vector<std::filesystem::path> fileName;
        fileName.push_back(test::TestFileManager::getInputFilePath("testINVEXP.xvg"));
        fileName.push_back(test::TestFileManager::getInputFilePath("testPRES.xvg"));
        fileName.push_back(test::TestFileManager::getInputFilePath("testINVEXP79.xvg"));
        fileName.push_back(test::TestFileManager::getInputFilePath("testERF.xvg"));
        fileName.push_back(test::TestFileManager::getInputFilePath("testERREST.xvg"));
        for (const auto& aName : fileName)
        {
            int        nrColumns;
            ExpfitData ed;
            ed.nrLines_   = read_xvg(aName, &tempValues, &nrColumns);
            ed.dt_        = tempValues[0][1] - tempValues[0][0];
            ed.startTime_ = tempValues[0][0];
            ed.endTime_   = tempValues[0][ed.nrLines_ - 1];
            for (int j = 0; j < ed.nrLines_; j++)
            {
                ed.x_.push_back(static_cast<real>(tempValues[0][j]));
                ed.y_.push_back(static_cast<real>(tempValues[1][j]));
            }
            data_.push_back(ed);

            // Free memory that was allocated in read_xvg
            for (int j = 0; j < nrColumns; j++)
            {
                sfree(tempValues[j]);
                tempValues[j] = nullptr;
            }
            sfree(tempValues);
            tempValues = nullptr;
        }
    }

    static void TearDownTestSuite() {}

    void test(int type, double result[], double tolerance, unsigned int testType)
    {
        int               nfitparm = effnNparams(type);
        gmx_output_env_t* oenv;

        if (testType >= data_.size())
        {
            GMX_THROW(InvalidInputError("testType out of range"));
        }
        output_env_init_default(&oenv);
        do_lmfit(data_[testType].nrLines_,
                 &(data_[testType].y_[0]),
                 nullptr,
                 data_[testType].dt_,
                 &(data_[testType].x_[0]),
                 data_[testType].startTime_,
                 data_[testType].endTime_,
                 oenv,
                 false,
                 type,
                 result,
                 0,
                 nullptr);
        output_env_done(oenv);
        checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, tolerance));
        checker_.checkSequenceArray(nfitparm, result, "result");
    }
};


// static var
std::vector<ExpfitData> ExpfitTest::data_;

// TODO calling test() leads to a fatal error, which we could in
// principle test for.

TEST_F(ExpfitTest, EffnEXP1)
{
    double param[] = { 25 };
    test(effnEXP1, param, 1e-5, 0);
}

TEST_F(ExpfitTest, EffnEXP2)
{
    double param[] = { 35, 0.5 };
    test(effnEXP2, param, 3e-5, 0);
}

TEST_F(ExpfitTest, EffnEXPEXP)
{
    double param[] = { 5, 0.5, 45 };
    test(effnEXPEXP, param, 1e-2, 0);
}

TEST_F(ExpfitTest, EffnEXP5)
{
    double param[] = { 0.5, 5, 0.5, 50, 0.002 };
    test(effnEXP5, param, 1e-2, 2);
}

TEST_F(ExpfitTest, EffnEXP7)
{
    double param[] = { 0.1, 2, 0.5, 30, 0.3, 50, -0.002 };
    test(effnEXP7, param, 1e-2, 2);
}

TEST_F(ExpfitTest, EffnEXP9)
{
    double param[] = { 0.4, 5, 0.2, 30, 0.1, 70, 0.2, 200, -0.05 };
    test(effnEXP9, param, 4e-2, 2);
}

TEST_F(ExpfitTest, EffnERF)
{
    double param[] = { 80, 120, 180, 5 };
    test(effnERF, param, 1e-1, 3);
}

TEST_F(ExpfitTest, EffnERREST)
{
    double param[] = { 1, 0.9, 100 };
    test(effnERREST, param, 5e-3, 4);
}

TEST_F(ExpfitTest, EffnVAC)
{
    double param[] = { 0.6, 0.1 };
    test(effnVAC, param, 0.05, 0);
}

TEST_F(ExpfitTest, EffnPRES)
{
    double param[] = { 0.6, 10, 7, 1, 0.25, 2 };
    test(effnPRES, param, 1e-4, 1);
}

#endif

} // namespace
} // namespace test
} // namespace gmx
