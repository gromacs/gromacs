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
 * Implements test of autocorrelation function routines
 *
 * \author Anders G&auml;rden&auml;s <anders.gardenas@gmail.com>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "gromacs/correlationfunctions/autocorr.h"

#include <cmath>

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/correlationfunctions/expfit.h"
#include "gromacs/fft/fft.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "correlationdataset.h"

namespace gmx
{
namespace test
{
namespace
{

//! Definition of pointer to class containing test data.
typedef std::unique_ptr<CorrelationDataSet> CorrelationDataSetPointer;

class AutocorrTest : public ::testing::Test
{
protected:
    static int                       nrFrames_;
    static CorrelationDataSetPointer data_;
    // Need raw pointer for passing this to C routines
    static t_pargs* tempArgs_;

    test::TestReferenceData    refData_;
    test::TestReferenceChecker checker_;

    // Use erefdataCreateMissing for creating new files
    AutocorrTest() : checker_(refData_.rootChecker())
    {
#if GMX_DOUBLE
        checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-6));
#else
        checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-3));
#endif
    }

    // Static initiation, only run once every test.
    static void SetUpTestSuite()
    {
        int         n        = 0;
        std::string fileName = "testCOS3.xvg";
        data_                = std::make_unique<CorrelationDataSet>(fileName);
        nrFrames_            = data_->getNrLines();
        tempArgs_            = add_acf_pargs(&n, nullptr);
    }

    static void TearDownTestSuite()
    {
        sfree(tempArgs_);
        tempArgs_ = nullptr;
        gmx_fft_cleanup();
    }

    void test(unsigned long mode, bool bNormalize)
    {
        bool              bAverage  = true;
        bool              bVerbose  = false;
        int               nrRestart = 1;
        int               dim       = getDim(mode);
        std::vector<real> result;

        for (int i = 0; i < nrFrames_; i++)
        {
            for (int m = 0; m < dim; m++)
            {
                result.push_back(data_->getValue(m, i));
            }
        }
        real* ptr = result.data();
        low_do_autocorr(nullptr,
                        nullptr,
                        nullptr,
                        nrFrames_,
                        1,
                        get_acfnout(),
                        &ptr,
                        data_->getDt(),
                        mode,
                        nrRestart,
                        bAverage,
                        bNormalize,
                        bVerbose,
                        data_->getStartTime(),
                        data_->getEndTime(),
                        effnNONE);

        double testResult = 0;
        for (int i = 0; i < get_acfnout(); i++)
        {
            testResult += result[i];
        }
        checker_.checkSequenceArray(get_acfnout(), ptr, "AutocorrelationFunction");
        checker_.checkReal(testResult, "Integral");
    }

    static int getDim(unsigned long type)
    {
        switch (type)
        {
            case eacNormal: return 1;
            case eacVector: return 3;
            case eacCos:
                return 1;
                // Several intended fall-throughs follow
            case eacRcross:
            case eacP0:
            case eacP1:
            case eacP2:
            case eacP3:
            case eacP4: return 3;
            case eacIden: return 1;
            default: GMX_RELEASE_ASSERT(false, "Invalid auto correlation option"); return -1;
        }
    }
};

int                       AutocorrTest::nrFrames_;
CorrelationDataSetPointer AutocorrTest::data_;
t_pargs*                  AutocorrTest::tempArgs_;

TEST_F(AutocorrTest, EacNormal)
{
    test(eacNormal, true);
}

TEST_F(AutocorrTest, EacNoNormalize)
{
    test(eacNormal, false);
}

TEST_F(AutocorrTest, EacCos)
{
    test(eacCos, true);
}

TEST_F(AutocorrTest, EacVector)
{
    test(eacVector, true);
}

TEST_F(AutocorrTest, EacRcross)
{
    test(eacRcross, true);
}

TEST_F(AutocorrTest, EacP0)
{
    test(eacP0, true);
}

TEST_F(AutocorrTest, EacP1)
{
    test(eacP1, true);
}

TEST_F(AutocorrTest, EacP2)
{
    test(eacP2, true);
}

TEST_F(AutocorrTest, EacP3)
{
    test(eacP3, true);
}

TEST_F(AutocorrTest, EacP4)
{
    test(eacP4, true);
}


} // namespace
} // namespace test
} // namespace gmx
