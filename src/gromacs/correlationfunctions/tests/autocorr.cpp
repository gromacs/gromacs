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
 * Implements test of autocorrelation function routines
 *
 * \author Anders G&auml;rden&auml;s <anders.gardenas@gmail.com>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "gromacs/correlationfunctions/autocorr.h"

#include <cmath>

#include <gtest/gtest.h>

#include "gromacs/correlationfunctions/expfit.h"
#include "gromacs/fft/fft.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/uniqueptr.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "correlationdataset.h"

namespace gmx
{
namespace
{

//! Definition of pointer to class containing test data.
typedef gmx_unique_ptr<CorrelationDataSet>::type CorrelationDataSetPointer;

class AutocorrTest : public ::testing::Test
{
    protected:

        static int                                  nrFrames_;
        static CorrelationDataSetPointer            data_;
        // Need raw pointer for passing this to C routines
        static t_pargs                            * tempArgs_;

        test::TestReferenceData                     refData_;
        test::TestReferenceChecker                  checker_;

        // Use erefdataCreateMissing for creating new files
        AutocorrTest( )
            : checker_(refData_.rootChecker())
        {
#ifdef GMX_DOUBLE
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-6));
#else
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-3));
#endif
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
            int         n        = 0;
            std::string fileName = "testCOS3.xvg";
            data_                = CorrelationDataSetPointer(new CorrelationDataSet(fileName));
            nrFrames_            = data_->getNrLines();
            tempArgs_            = add_acf_pargs(&n, NULL);
        }

        static void TearDownTestCase()
        {

            sfree(tempArgs_);
            tempArgs_ = NULL;
            gmx_fft_cleanup();
        }

        void test(unsigned long mode)
        {
            bool              bAverage      = false;
            bool              bNormalize    = true;
            bool              bVerbose      = false;
            int               nrRestart     = 1;
            int               dim           = getDim(mode);
            std::vector<real> result;

            for (int i = 0; i < nrFrames_; i++)
            {
                for (int m = 0; m < dim; m++)
                {
                    result.push_back(data_->getValue(m, i));
                }
            }
            real *ptr = static_cast<real*>(&(result[0]));
            low_do_autocorr(0, 0, 0,   nrFrames_, 1,
                            get_acfnout(), &ptr, data_->getDt(), mode,
                            nrRestart, bAverage, bNormalize,
                            bVerbose, data_->getStartTime(), data_->getEndTime(),
                            effnNONE);

            double testResult = 0;
            for (int i = 0; i < nrFrames_; i++)
            {
                testResult += result[i];
            }
            checker_.checkSequenceArray(nrFrames_, ptr,
                                        "AutocorrelationFunction");
            checker_.checkReal(testResult, "Integral");
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
                    GMX_RELEASE_ASSERT(false, "Invalid auto correlation option");
                    return -1;
            }

        }

};

int                         AutocorrTest::nrFrames_;
CorrelationDataSetPointer   AutocorrTest::data_;
t_pargs                   * AutocorrTest::tempArgs_;

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
