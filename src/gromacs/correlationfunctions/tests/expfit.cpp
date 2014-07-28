/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * Implements test of exponential fitting routines
 *
 * \author Anders G&auml;rden&auml;s <anders.gardenas@gmail.com>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "gromacs/correlationfunctions/expfit.h"

#include <cmath>

#include <gtest/gtest.h>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

//! Number of data files for testing.
#define expTestNrTypes 3

namespace gmx
{

namespace
{

class ExpfitTest : public ::testing::Test
{

    protected:
        static int                 nrLines_;
        static std::vector<real>   values_[expTestNrTypes];
        static int                 nrColumns_;
        static std::vector<real>   standardDev_;
        static real                startTime_;
        static real                endTime_;
        static real                timeDeriv_;
        test::TestReferenceData    refData_;
        test::TestReferenceChecker checker_;
        ExpfitTest( )
            : checker_(refData_.rootChecker())
        {
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
            double   ** tempValues;
            std::string fileName[expTestNrTypes];
            fileName[0] = test::TestFileManager::getInputFilePath("testINVEXP.xvg");
            fileName[1] = test::TestFileManager::getInputFilePath("testPRES.xvg");
            fileName[2] = test::TestFileManager::getInputFilePath("testEXP.xvg");
            for (int i = 0; i < expTestNrTypes; i++)
            {
                const char * name = fileName[i].c_str();
                // TODO: this assumes all files have the same length.
                nrLines_     = read_xvg(name, &tempValues, &nrColumns_);

                // Generating standard deviation
                if (i == 0)
                {
                    double fac = 1.0/nrLines_;
                    for (int j = 0; j < nrLines_; j++)
                    {
                        standardDev_.push_back(fac);
                    }
                    timeDeriv_ = tempValues[0][1] - tempValues[0][0];
                    startTime_ = tempValues[0][0];
                    endTime_   = tempValues[0][nrLines_-1];
                }

                for (int j = 0; j  < nrLines_; j++)
                {
                    values_[i].push_back((real)tempValues[1][j]);
                }

                // Free memory that was allocated in read_xvg
                for (int i = 0; i < nrColumns_; i++)
                {
                    sfree(tempValues[i]);
                    tempValues[i] = NULL;
                }
                sfree(tempValues);
                tempValues = NULL;
            }
        }

        static void TearDownTestCase()
        {
        }

        void test(int type, double result[], double tolerance, int testType)
        {
            int     nfitparm = effnNparams(type);

            do_lmfit(nrLines_, &values_[testType][0], &standardDev_[0], timeDeriv_,
                     NULL, startTime_, endTime_, NULL, false, type, result, 0);

            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, tolerance));
            checker_.checkSequenceArray(nfitparm, result, "result");
        }
};


//static var
int               ExpfitTest::nrLines_;
//cppcheck-suppress arrayIndexOutOfBounds fixed in 1.68-dev
std::vector<real> ExpfitTest::values_[expTestNrTypes];
int               ExpfitTest::nrColumns_;
std::vector<real> ExpfitTest::standardDev_;
real              ExpfitTest::startTime_;
real              ExpfitTest::endTime_;
real              ExpfitTest::timeDeriv_;

TEST_F (ExpfitTest, EffnEXP1) {
    double  param[] = {25};
    test(effnEXP1, param, 1e-6, 0);
}

TEST_F (ExpfitTest, EffnEXP2) {
    double  param[] = {35, 0.5};
    test(effnEXP2, param, 1e-6, 0);
}

TEST_F (ExpfitTest, EffnEXP3) {
    double param[] = {45, 0.5, 5};
    test(effnEXP3, param, 1e-4, 0);
}

TEST_F (ExpfitTest, EffnEXP5) {
    double  param[] = {0.5, 5, 0.5, 50, 0.002};
    test(effnEXP5, param, 1e-4, 0);
}

TEST_F (ExpfitTest, EffnEXP7) {
    double  param[] = {0.5, 5, -0.02, 0.5, 0.5, 50, -0.002};
    test(effnEXP7, param, 1e-4, 0);
}

TEST_F (ExpfitTest, EffnEXP9) {
    double  param[] = {2, 1200, -1, 300, 0.7, 70, 0.5, 6, -0.5};
    test(effnEXP9, param, 4e-2, 0);
}

TEST_F (ExpfitTest, EffnERF) {
    double  param[] = {0.5, 0.5, 0.5, 1};
    test(effnERF, param, 1e-2, 0);
}

TEST_F (ExpfitTest, EffnERREST) {
    double  param[] = {0.5, 0.7, 0.3};
    test(effnERREST, param, 1e-4, 2);
}

TEST_F (ExpfitTest, EffnVAC) {
    double param[] = {0.5, 0.05};
    test(effnVAC, param, 1e-4, 0);
}

TEST_F (ExpfitTest, DISABLED_EffnPRES) {
    //TODO: This test is prodocues NaNs and INFs. Fix and then reactivate.
    double param[] = {0, 10, 4, 1, 0.5, 1};
    test(effnPRES, param, 1e-4, 1);
}

}

}
