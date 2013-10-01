/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
 * Implements test routines from scattering factor IO
 *
 * \author Daccid van der Spoel
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include <gtest/gtest.h>
#include <cstring>

#include "gromacs/waxsdebye/scattering_factors.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
//#include "testutils/testfilemanager.h"

namespace gmx
{
namespace
{
class ScatteringFactorTableTest : public ::testing::Test
{

    protected:
        static gmx::ScatteringFactorTable sft;

        //static init, only run once every test.
        static void SetUpTestCase(const char *datafile)
        {
            sft.read(datafile);
        }

        test::TestReferenceData    refData;
        test::TestReferenceChecker checker;


        ScatteringFactorTableTest( )
            : refData(test::erefdataCreateMissing), checker(refData.rootChecker())
        {

#ifdef GMX_DOUBLE
            checker.setDefaultTolerance(test::relativeRealTolerance(1, 20000));
#else
            checker.setDefaultTolerance(test::relativeRealTolerance(1, 5000));
#endif
        }

        static void TearDownTestCase()
        {
        }

        void test(const char *datafile)
        {
            SetUpTestCase(datafile);

            //int     testType = getTestType(type);
            /*
               real  * values   = new real[data[0]->getNrLines()];
               for (int i = 0; i < data[0]->getNrLines(); i++)
               {
                values[i] = data[testType]->getValue(i);

                }*/
            //init data and cheker might be a beter plays to init them

            //testResult(result, nfitparm);
            /* delete [] values; */
        }

        void testResult(real result[], int nrFitParam)
        {

            std::string testName = "result";
            checker.checkSequenceArray(nrFitParam, result, testName.c_str());

        }

        int getTestType(unsigned long type)
        {
            return (int)type;
        }

        std::string convertInt(int number)
        {
            std::stringstream ss; //create a stringstream
            ss << number;         //add number to the stream
            return ss.str();      //return a string with the contents of the stream
        }

};

TEST_F(ScatteringFactorTableTest, MartiniDS)
{
    const char *sfactorFile = "sfactor_martini_ds_Fourier.xml";
    test(sfactorFile);
}

TEST_F(ScatteringFactorTableTest, FileIO)
{

    //    sft.write("sfactor_martini_ds.out");
}

}

}
