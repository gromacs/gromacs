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
#include "gromacs/fft/fft.h"
#include <testutils/testfilemanager.h>

#include "gromacs/utility/smalloc.h"
#include "gromacs/correlationfunctions/autocorr.h"


class Autocorr : public ::testing::Test
{
    protected:
        static std::string fileName;
        static int         nrFrames;
        static real      * result;
        static real      * values3D;
        static int         nrColums;
        static real        startTime;
        static real        endTime;
        static real        dt;
        static std::string resultName;
        static t_pargs   * tempArgs;

        Autocorr( )
        {
        }

        //static init, only runed onecs every test.
        static void SetUpTestCase()
        {
            int       n = 0;
            double ** values;
            int       nrLines;
            tempArgs    = add_acf_pargs(&n, NULL);
            fileName    = gmx::test::TestFileManager::getInputFilePath("testCOS3.xvg");

            nrLines     = read_xvg(fileName.c_str(), &values, &nrColums );
            nrFrames    = nrLines/3;
            values3D    = new real[nrLines];
            result      = new real[nrLines];




            //Comverting double to real
            for (int i = 0; i  < nrLines; i++)
            {
                values3D[i]   = (real)values[1][i];

            }
            dt        =  values[0][1] - values[0][0];
            startTime = values[0][0];
            endTime   = values[0][nrLines-1];




            //alocated in read_xvg
            for (int i = 0; i < nrColums; i++)
            {
                sfree(values[i]);
                values[i] = NULL;
            }
            sfree(values);
            values = NULL;


        }

        void SetUp( )
        {
            //Seting up ref test
        }

        static void TearDownTestCase()
        {
            delete [] values3D;
            values3D = NULL;
            delete [] result;
            result = NULL;
            sfree(tempArgs);
            tempArgs = NULL;

            gmx_fft_cleanup();
        }

        void TearDown( )
        {
        }


        void test(unsigned long mode)
        {
            bool startSame     = true;
            bool normalize     = false;
            bool printConsole  = true;    int nrRestart = 1;
            int  nrFitingParam = 0; // get_acffitfn()
            int  dim           = getDim(mode);
            if (dim == 1)
            {
                for (int i = 0; i < nrFrames; i++)
                {
                    result[i] = values3D[i*3];
                }
            }
            else
            {
                for (int i = 0; i < nrFrames*dim; i++)
                {
                    result[i] = values3D[i];
                }
            }
            low_do_autocorr(0, 0,
                            0,   nrFrames, 1,
                            get_acfnout(), &result, dt, mode,
                            nrRestart, normalize, startSame,
                            printConsole, startTime, endTime,
                            nrFitingParam);

            gmx::test::TestReferenceData    data(gmx::test::erefdataCreateMissing);
            gmx::test::TestReferenceChecker checker(data.rootChecker());
            testResult(result, checker, "result");

        }





        void testResult(real result[], gmx::test::TestReferenceChecker checker, std::string testName)
        {
            real testResult = 0;
            for (int i = 0; i < nrFrames; i++)
            {
                testResult +=  result[i];
            }
            double percentBound = 0.01;
            checker.setDefaultTolerance(gmx::test::FloatingPointTolerance(std::abs(testResult) * percentBound, -1, false));
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
                    printf("invalid Autocorrelation type");
                    return 1;
            }
        }

};


//static var
std::string Autocorr::fileName;
int         Autocorr::nrFrames;

real      * Autocorr::result;
real      * Autocorr::values3D;
int         Autocorr::nrColums;
real        Autocorr::startTime;
real        Autocorr::endTime;
real        Autocorr::dt;
t_pargs   * Autocorr::tempArgs;



TEST_F (Autocorr, EacNormal) {
    test(eacNormal);
}


TEST_F (Autocorr, EacCos) {
    test(eacCos);
}

TEST_F (Autocorr, EacVector) {
    test(eacVector);
}

TEST_F (Autocorr, EacRcross) {
    test(eacRcross);
}

TEST_F (Autocorr, EacP0) {
    test(eacP0);
}

TEST_F (Autocorr, EacP1) {
    test(eacP1);
}

TEST_F (Autocorr, EacP2) {
    test(eacP2);
}

TEST_F (Autocorr, EacP3) {
    test(eacP3);
}

TEST_F (Autocorr, EacP4) {
    test(eacP4);
}


TEST_F (Autocorr, EacIden) {
    //test(eacIden);
    //not supported
}
