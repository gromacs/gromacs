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


#define expTestNrTypes 3

#include <cmath>
#include <oenv.h>
#include <gromacs/fileio/xvgr.h>
#include <gtest/gtest.h>
#include <sstream>
#include <testutils/refdata.h>
#include <testutils/testasserts.h>
#include <testutils/testfilemanager.h>

#include "gromacs/utility/smalloc.h"
#include "gromacs/correlationfunctions/expfit.h"

class Expfit : public ::testing::Test
{

    protected:
        static int         nrLines;
        static real      * values[expTestNrTypes];
        static int         nrColums;
        static real      * standardDiv;
        static real        startTime;
        static real        endTime;
        static real        timeDirev;

        //static init, only runed onecs every test.
        static void SetUpTestCase()
        {
            double   ** tempValues;
            std::string fileName[expTestNrTypes];
            fileName[0] = gmx::test::TestFileManager::getInputFilePath("testINVEXP.xvg");
            fileName[1] = gmx::test::TestFileManager::getInputFilePath("testPRES.xvg");
            fileName[2] = gmx::test::TestFileManager::getInputFilePath("testEXP.xvg");
            for (int i = 0; i < expTestNrTypes; i++)
            {

                const char * name = fileName[i].c_str();
                nrLines     = read_xvg(name, &tempValues, &nrColums);

                //generating standard div
                if (i == 0)
                {
                    real fac = 1.0/((real)nrLines);
                    standardDiv = new real[nrLines];
                    for (int j = 0; j < nrLines; j++)
                    {
                        standardDiv[j] = fac;
                    }
                    timeDirev =  tempValues[0][1] - tempValues[0][0];
                    startTime = tempValues[0][0];
                    endTime   = tempValues[0][nrLines-1];
                }


                values[i] = new real[nrLines];
                for (int j = 0; j  < nrLines; j++)
                {
                    values[i][j]   = (real)tempValues[1][j];
                }

                //alocated in read_xvg
                for (int i = 0; i < nrColums; i++)
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
            delete [] standardDiv;
            standardDiv = NULL;
            for (int i = 0; i < expTestNrTypes; i++)
            {
                delete [] values[i];
                values[i] = NULL;
            }
        }



        void test(int type, real result[])
        {
            int     nfitparm = effnNparams(type);
            int     testType = getTestType(type);

            //init data and cheker might be a beter plays to init them
            gmx::test::TestReferenceData    data(gmx::test::erefdataCreateMissing);
            gmx::test::TestReferenceChecker checker(data.rootChecker());
            do_lmfit(nrLines, values[testType], standardDiv, timeDirev, NULL, startTime, endTime, NULL, false, type, result, 0);
            testResult(result, nfitparm, checker, "Result");
        }




        void testResult(real result[], int nrFitParam, gmx::test::TestReferenceChecker checker, std::string testName)
        {
            double percentBound = 0.01;

            for (int i = 0; i < nrFitParam; i++)
            {
                checker.setDefaultTolerance(gmx::test::FloatingPointTolerance(std::abs(result[i]) * percentBound, -1, false));
                checker.checkReal(result[i], (testName + convertInt(i)).c_str());
            }
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
                    printf("invalid expfit type \n");
                    return 0;
            }
        }

        //std::to_string does not compile
        std::string convertInt(int number)
        {
            std::stringstream ss; //create a stringstream
            ss << number;         //add number to the stream
            return ss.str();      //return a string with the contents of the stream
        }

};


//static var
int         Expfit::nrLines;
real      * Expfit::values[expTestNrTypes];
int         Expfit::nrColums;
real      * Expfit::standardDiv;
real        Expfit::startTime;
real        Expfit::endTime;
real        Expfit::timeDirev;

TEST_F (Expfit, EffnEXP1) {
    real  param[] = {25};
    test(effnEXP1, param);
}

TEST_F (Expfit, EffnEXP2) {
    real  param[] = {35, 0.5};
    test(effnEXP2, param);
}

TEST_F (Expfit, EffnEXP3) {
    real param[] = {45, 0.5, 5};
    test(effnEXP3, param);
}

TEST_F (Expfit, EffnEXP5) {
    real  param[] = {0.5, 5, 0.5, 50, 0.002};
    test(effnEXP5, param);
}

TEST_F (Expfit, EffnEXP7) {
    real  param[] = {0.5, 5, -0.02, 0.5, 0.5, 50, -0.002};
    test(effnEXP7, param);
}

TEST_F (Expfit, EffnEXP9) {
    real  param[] = {0.18, 800, -0.2, 161, 0.7, 60, 0.5, 5, 0.1};
    test(effnEXP9, param);
}

TEST_F (Expfit, EffnERF) {
    real  param[] = {0.5, 0.5, 0.5, 1};
    test(effnERF, param);
}

TEST_F (Expfit, EffnERREST) {
    real  param[] = {0.5, 0.7, 0.3};
    test(effnERREST, param);
}

TEST_F (Expfit, EffnVAC) {
    real param[] = {0.5, 0.05};
    test(effnVAC, param);
}

TEST_F (Expfit, EffnPRES) {
    real param[] = {0, 10, 4, 1, 0.5, 1};
    test(effnPRES, param);
}
