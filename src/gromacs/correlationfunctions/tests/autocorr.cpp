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

#include <cmath>
#include <sstream>
#include "gtest/gtest.h"
#include "gromacs/utility/gmxassert.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "gromacs/fft/fft.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "correlationDataSet.h"

namespace gmx
{
namespace
{
class AutocorrTest : public ::testing::Test
{
    protected:

        static int                    nrFrames;
        static   CorrelationDataSet * data;
        static  t_pargs             * tempArgs;




        gmx::test::TestReferenceData    refData;
        gmx::test::TestReferenceChecker checker;


        AutocorrTest( )
            : refData(gmx::test::erefdataCompare), checker(refData.rootChecker())
        {

#ifdef GMX_DOUBLE
            checker.setDefaultTolerance(gmx::test::relativeRealTolerance(1, 20000));
#else
            checker.setDefaultTolerance(gmx::test::relativeRealTolerance(1, 5000));
#endif
        }

        //static init, only runed onecs every test.
        static void SetUpTestCase()
        {
            int         n        = 0;
            std::string fileName = "testCOS3.xvg";
            data        = new CorrelationDataSet(fileName);
            nrFrames    = data->getNrLines()/3;
            tempArgs    = add_acf_pargs(&n, NULL);

        }


        static void TearDownTestCase()
        {

            sfree(tempArgs);
            tempArgs = NULL;
            delete data;
            gmx_fft_cleanup();
        }


        void test(unsigned long mode)
        {
            bool   startSame     = true;
            bool   normalize     = false;
            bool   printConsole  = true;    int nrRestart = 100;
            int    nrFitingParam = 0; // get_acffitfn()
            int    dim           = getDim(mode);
            real * result        = new real[nrFrames*dim];
            if (dim == 1)
            {
                for (int i = 0; i < nrFrames; i++)
                {
                    result[i] = data->getValue(i*3);
                }
            }
            else
            {
                for (int i = 0; i < nrFrames*dim; i++)
                {
                    result[i] = data->getValue(i);
                }
            }
            low_do_autocorr(0, 0, 0,   nrFrames, 1,
                            get_acfnout(), &result, data->getDt(), mode,
                            nrRestart, normalize, startSame,
                            printConsole, data->getStartTime(), data->getEndTime(),
                            nrFitingParam);


            testResult(result);
            delete [] result;
        }

        void testResult(real result[])
        {
            real testResult = 0;
            for (int i = 0; i < nrFrames; i++)
            {
                testResult +=  result[i];
            }
            std::string testName = "IntegralSize";
            checker.checkReal(testResult, testName.c_str());
        }

        int getDim(unsigned long type)
        {
            switch (type)
            {
                case eacNormal:
                    return 1;
                case eacVector:
                    return 3;
                case eacCos:
                    return 1;
                case eacRcross:
                    return 3;
                case eacP0:
                    return 3;
                case eacP1:
                    return 3;
                case eacP2:
                    return 3;
                case eacP3:
                    return 3;
                case eacP4:
                    return 3;
                case eacIden:
                    return 1;
                default:
                    GMX_RELEASE_ASSERT(true, "Invalid auto correlation option");
                    return -1;
            }
// local variables


        }

};

int                  AutocorrTest::nrFrames;
CorrelationDataSet * AutocorrTest::data;
t_pargs            * AutocorrTest::tempArgs;

TEST_F (AutocorrTest, EacNormal)
{
    test(eacNormal);
}

TEST_F (AutocorrTest, EacCos)
{
    test(eacCos);
}

TEST_F (AutocorrTest, EacVector)
{
    test(eacVector);
}

TEST_F (AutocorrTest, EacRcross)
{
    test(eacRcross);
}

TEST_F (AutocorrTest, EacP0)
{
    test(eacP0);
}

TEST_F (AutocorrTest, EacP1)
{
    test(eacP1);
}

TEST_F (AutocorrTest, EacP2)
{
    test(eacP2);
}

TEST_F (AutocorrTest, EacP3)
{
    test(eacP3);
}

TEST_F (AutocorrTest, EacP4)
{
    test(eacP4);
}


}

}
