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
#include <oenv.h>
#include <gromacs/fileio/xvgr.h>
#include <gtest/gtest.h>
#include <sstream>
#include <testutils/refdata.h>
#include <testutils/testasserts.h>

#include "gromacs/utility/smalloc.h"
#include "../expfit.h"


class Expfit : public ::testing::Test
{
    protected:
        static std::string fileName;
        static int         nrLines;
        static real      * values[2];
        static int         nrColums;
        static real      * standardDiv;
        static real        startTime;
        static real        endTime;
        static real        timeDirev;


        output_env_t oenv;

        Expfit( )
        {
            oenv = NULL;
            //
            // initialization code here
        }

        //static init, only runed onecs every test.
        static void SetUpTestCase()
        {
            double ** values2;
            fileName    = gmx::test::getReferenceDataPath() + "/test1.xvg";
            nrColums    = 2;
            const char * name = fileName.c_str();
            nrLines     = read_xvg(name, &values2, &nrColums );
            standardDiv = new real[nrLines];

            values[0] = new real[nrLines];
            values[1] = new real[nrLines];


            //generating standard div
            real fac = 1.0/((real)nrLines);

            //Comverting double to real
            for (int i = 0; i  < nrLines; i++)
            {
                values[0][i]   = (real)values2[0][i];
                values[1][i]   = (real)values2[1][i];
                standardDiv[i] = fac;
            }

            timeDirev =  values[0][1] - values[0][0];
            startTime = values[0][0];
            endTime   = values[0][nrLines-1];

            //alocated in read_xvg
            for (int i = 0; i < nrColums; i++)
            {
                sfree(values2[i]);
                values2[i] = NULL;
            }
            sfree(values2);
            values2 = NULL;

        }

        void SetUp( )
        {
            //Seting up ref test
            //  output_env_init_default(&oenv);



            // code here will execute just before the test ensues
        }

        static void TearDownTestCase()
        {
            delete [] standardDiv;
            standardDiv = NULL;
            delete [] values[0];
            values[0] = NULL;
            delete [] values[1];
            values[1] = NULL;
        }

        void TearDown( )
        {
        }

        ~Expfit( )
        {
        }

        void test(int Type, real result[])
        {
            int     nfitparm = effnNparams(Type);
            /*    real  * result   = new real[nfitparm];
                for (int i = 0; i < nfitparm; i++)
                {
                    result[i] = 0.5;
               }*/
            gmx::test::TestReferenceData    data(gmx::test::erefdataCreateMissing);
            gmx::test::TestReferenceChecker checker(data.rootChecker());
            testResult(result, nfitparm, checker, "Param");
            do_lmfit(nrLines, values[1], standardDiv, timeDirev, values[0], startTime, endTime, oenv,
                     false, Type, result, 0);
            testResult(result, nfitparm, checker, "Result");
            //delete [] result;
        }


        void testResult(real result[], int size, gmx::test::TestReferenceChecker checker, std::string testName)
        {
            double percentBound = 0.05;
            for (int i = 0; i < size; i++)
            {

                checker.setDefaultTolerance(gmx::test::FloatingPointTolerance(
                                                    abs(result[i]) * percentBound, 4, false));
                checker.checkReal(result[i], (testName + convertInt(i)).c_str());
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
std::string Expfit::fileName;
int         Expfit::nrLines;
real      * Expfit::values[2];
int         Expfit::nrColums;
real      * Expfit::standardDiv;
real        Expfit::startTime;
real        Expfit::endTime;
real        Expfit::timeDirev;

TEST_F (Expfit, effnEXP1) {
    real  param[] = {25};
    test(effnEXP1, param);
}

TEST_F (Expfit, effnEXP2) {
    real  param[] = {35, 0.5};
    test(effnEXP2, param);
}

TEST_F (Expfit, effnEXP3) {
    real param[] = {45, 0.5, 5};
    test(effnEXP3, param);
}

TEST_F (Expfit, effnEXP5) {
    real  param[] = {0.5, 5, 0.5, 50, 0.002};
    test(effnEXP5, param);
}

TEST_F (Expfit, effnEXP7) {
    real  param[] = {0.5, 5, -0.02, 0.5, 0.5, 50, -0.002};
    test(effnEXP7, param);
}

TEST_F (Expfit, effnEXP9) {
    real  param[] = {0.18, 800, -0.2, 161, 0.7, 60, 0.5, 5, 0.1};
    test(effnEXP9, param);
}

TEST_F (Expfit, effnERF) {
    real  param[] = {-54, 54, 0.5, 0.5};
    test(effnERF, param);
}

TEST_F (Expfit, effnERREST) {
    real  param[] = {0.5, 0.5, 0.5};
    test(effnERREST, param);
}

TEST_F (Expfit, effnVAC) {
    real param[] = {0.5, 0.05};
    test(effnVAC, param);
}
