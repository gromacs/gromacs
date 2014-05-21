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
#include <oenv.h>
#include <gromacs/fileio/xvgr.h>
#include <gtest/gtest.h>

#include "../expfit.h"
#include <testutils/refdata.h>


class Expfit : public ::testing::Test
{
    protected:
        static std::string fileName;
        static int         nrLines;
        static float    ** values;
        static int         nrColums;
        static float     * standardDiv;
        static float       startTime;
        static float       endTime;
        static float       timeDirev;


        output_env_t oenv;

        Expfit( )
        {

            // initialization code here
        }


        static void SetUpTestCase()
        {
            fileName = gmx::test::getReferenceDataPath() + "/test1.xvg";
            nrColums = 2;
            double ** values2;
            nrLines     = read_xvg(fileName.c_str(), &values2, &nrColums );
            standardDiv = new float[nrLines];

            values    = new float*[2];
            values[0] = new float[nrLines];
            values[1] = new float[nrLines];


            //generating standard div
            float fac = 1.0/((real)nrLines);


            for (int i = 0; i  < nrLines; i++)
            {
                values[0][i]   = (float)values2[0][i];
                values[1][i]   = (float)values2[1][i];
                standardDiv[i] = fac;
            }

            timeDirev =  values[0][1] - values[0][0];
            startTime = values[0][0];
            endTime   = values[0][nrLines-1];
        }

        void SetUp( )
        {
            //Seting up ref test
            //  output_env_init_default(&oenv);



            // code here will execute just before the test ensues
        }

        static void TearDownTestCase()
        {
            delete standardDiv;
            standardDiv = NULL;
            delete values[0];
            values[0] = NULL;
            delete values[1];
            values[1] = NULL;
            delete values;
            values = NULL;
        }

        void TearDown( )
        {
        }

        ~Expfit( )
        {
// cleanup any pending stuff, but no exceptions allowed
        }

        // put in any custom data members that you need


        void test(int Type)
        {
            int     nfitparm = effnNparams(Type);
            float * result   = new float[nfitparm];
            for (int i = 0; i < nfitparm; i++)
            {
                result[i] = 1;
            }
            gmx::test::TestReferenceData    data(gmx::test::erefdataUpdateAll);
            gmx::test::TestReferenceChecker checker(data.rootChecker());
            testParam(result, nfitparm, checker);
            do_lmfit(nrLines, values[1], standardDiv, timeDirev, values[0], startTime, endTime, oenv,
                     false, Type, result, 0);
            testResult(result, nfitparm, checker);

        }


        void testParam(float params[], int size, gmx::test::TestReferenceChecker checker)
        {
            for (int i = 0; i < size; i++)
            {
                checker.checkFloat(params[i], ("Param"+ std::to_string(i)).c_str());
            }
        }

        void testResult(float result[], int size, gmx::test::TestReferenceChecker checker)
        {
            for (int i = 0; i < size; i++)
            {
                checker.checkFloat(result[i], ("Result" + std::to_string(i)).c_str());
            }
        }

};


//static var
std::string Expfit::fileName;
int         Expfit::nrLines;
float    ** Expfit::values;
int         Expfit::nrColums;
float     * Expfit::standardDiv;
float       Expfit::startTime;
float       Expfit::endTime;
float       Expfit::timeDirev;

TEST_F (Expfit, effnEXP1) {
    test(effnEXP1);
}

TEST_F (Expfit, effnEXP2) {
    test(effnEXP2);
}

TEST_F (Expfit, effnEXP3) {
    test(effnEXP3);
}

TEST_F (Expfit, effnEXP5) {
    test(effnEXP5);
}

TEST_F (Expfit, effnEXP7) {
    test(effnEXP7);
}

TEST_F (Expfit, effnEXP9) {
    test(effnEXP9);
}

TEST_F (Expfit, effnERF) {
    test(effnERF);
}

TEST_F (Expfit, effnERREST) {
    test(effnERREST);
}

TEST_F (Expfit, effnVAC) {
    test(effnVAC);
}
