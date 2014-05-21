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
#include "gromacs/utility/gmxassert.h"
#include "gtest/gtest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/correlationfunctions/expfit.h"
#include "correlationDataSet.h"


namespace gmx
{
namespace
{

#define expTestNrTypes 3

class ExpfitTest : public ::testing::Test
{

    protected:
        static CorrelationDataSet * data[3];
        static real              *  standardDiv;

        //static init, only runed onecs every test.
        static void SetUpTestCase()
        {

            std::string fileName = "testINVEXP.xvg";
            data[0]  = new CorrelationDataSet(fileName);
            fileName = "testPRES.xvg";
            data[1]  = new CorrelationDataSet(fileName);
            fileName = "testEXP.xvg";
            data[2]  = new CorrelationDataSet(fileName);


            real fac = 1.0/((real)data[0]->getNrLines());
            standardDiv = new real[data[0]->getNrLines()];
            for (int j = 0; j < data[0]->getNrLines(); j++)
            {
                standardDiv[j] = fac;
            }

        }



        gmx::test::TestReferenceData    refData;
        gmx::test::TestReferenceChecker checker;


        ExpfitTest( )
            : refData(gmx::test::erefdataCompare), checker(refData.rootChecker())
        {

#ifdef GMX_DOUBLE
            checker.setDefaultTolerance(gmx::test::relativeRealTolerance(1, 40000));
#else
            checker.setDefaultTolerance(gmx::test::relativeRealTolerance(1, 20000));
#endif
        }




        static void TearDownTestCase()
        {
            delete [] standardDiv;
            standardDiv = NULL;
            for (int i = 0; i < 3; i++)
            {
                delete  data[i];

            }
        }

        void test(int type, real result[])
        {
            int     nfitparm = effnNparams(type);
            int     testType = getTestType(type);
            real  * values   = new real[data[0]->getNrLines()];
            for (int i = 0; i < data[0]->getNrLines(); i++)
            {
                values[i] = data[testType]->getValue(i);

            }
            //init data and cheker might be a beter plays to init them


            do_lmfit(data[0]->getNrLines(), values, standardDiv, data[0]->getDt(), NULL,
                     data[0]->getStartTime(), data[0]->getEndTime(), NULL, false, type, result, 0);


            testResult(result, nfitparm);
            delete [] values;
        }

        void testResult(real result[], int nrFitParam)
        {

            std::string testName = "result";
            checker.checkSequenceArray(nrFitParam, result, testName.c_str());

        }





        int getTestType(unsigned long type)
        {
            /*
               type 0 exp test
               type 1 PRES test
               type 2 ERF test
             */
            switch (type)
            {
                case effnEXP1:
                    return 0;
                case effnEXP2:
                    return 0;
                case effnEXP3:
                    return 0;
                case effnEXP5:
                    return 0;
                case effnEXP7:
                    return 0;
                case effnEXP9:
                    return 0;
                case effnVAC:
                    return 0;
                case effnERF:
                    return 0;
                case effnERREST:
                    return 2;
                case effnPRES:
                    return 1;
                default:
                    GMX_RELEASE_ASSERT(true, "Invalid expfit option");
                    return -1;
            }
        }


        std::string convertInt(int number)
        {
            std::stringstream ss; //create a stringstream
            ss << number;         //add number to the stream
            return ss.str();      //return a string with the contents of the stream
        }

};


CorrelationDataSet * ExpfitTest::data[3];
real               * ExpfitTest::standardDiv;






TEST_F (ExpfitTest, EffnEXP1)
{
    real  param[] = {25};
    test(effnEXP1, param);
}

TEST_F (ExpfitTest, EffnEXP2)
{
    real  param[] = {35, 0.5};
    test(effnEXP2, param);
}

TEST_F (ExpfitTest, EffnEXP3)
{
    real param[] = {45, 0.5, 5};
    test(effnEXP3, param);
}

TEST_F (ExpfitTest, EffnEXP5)
{
    real  param[] = {0.5, 5, 0.5, 50, 0.002};
    test(effnEXP5, param);
}

TEST_F (ExpfitTest, EffnEXP7)
{
    real  param[] = {0.5, 5, -0.02, 0.5, 0.5, 50, -0.002};
    test(effnEXP7, param);
}

TEST_F (ExpfitTest, EffnEXP9)
{
    real  param[] = {0.18, 800, -0.2, 161, 0.7, 60, 0.5, 5, 0.1};
    test(effnEXP9, param);
}

TEST_F (ExpfitTest, EffnERF)
{
    real  param[] = {0.5, 0.5, 0.5, 1};
    test(effnERF, param);
}

TEST_F (ExpfitTest, EffnERREST) {
    real  param[] = {0.5, 0.7, 0.3};
    test(effnERREST, param);
}

TEST_F (ExpfitTest, EffnVAC)
{
    real param[] = {0.5, 0.05};
    test(effnVAC, param);
}

TEST_F (ExpfitTest, EffnPRES)
{
    real param[] = {0, 10, 4, 1, 0.5, 1};
    test(effnPRES, param);
}

}

}
