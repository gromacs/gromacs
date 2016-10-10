/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Implements PME B-spline moduli computation tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include <gmock/gmock.h>

#include "gromacs/ewald/pme-internal.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "pmeabstracttest.h"

/*! \brief Moduli algorithm type */
enum class ModuliType
{
    PME,
    P3M
};
// When adding new members to enums, don't forget to alter the test instantiations!

/*! \brief Convenience typedef of input parameters - grid dimensions, PME interpolation order, moduli type */
typedef std::tuple<gmx::IVec, int, ModuliType> BSplineModuliInputParameters;
/*! \brief Convenience typedef of test parameters - input parameters + tester function pointer */
typedef std::tuple<BSplineModuliInputParameters, TesterFuncPtr> BSplineModuliTestParameters;

/*! \brief Test fixture for testing PME B-spline moduli creation */
class PmeBSplineModuliTest : public ::testing::TestWithParam<BSplineModuliTestParameters>, public AbstractPmeTest
{

    //! PME spline interpolation order
    int         pmeOrder_;
    //! PME grid size
    gmx::IVec   gridSize_;
    //! Moduli algorithm type - see above
    ModuliType  moduliType_;
    public:
        PmeBSplineModuliTest(){}
        /*! \brief An initializer */
        void SetUp()
        {
            /* Getting the input */
            BSplineModuliInputParameters parameters;
            std::tie(parameters, testerFunction_)       = GetParam();
            std::tie(gridSize_, pmeOrder_, moduliType_) = parameters;

            /* Storing the input where it's needed */
            inputRec_->nkx         = gridSize_[XX];
            inputRec_->nky         = gridSize_[YY];
            inputRec_->nkz         = gridSize_[ZZ];
            inputRec_->coulombtype = (moduliType_ == ModuliType::P3M) ? eelP3M_AD : eelPME;
            inputRec_->pme_order   = pmeOrder_;

            /* Describing the test uniquely */
            std::string testDescription = gmx::formatString("Testing B-spline moduli creation (%s) for PME order %d, grid size %d %d %d",
                                                            (moduliType_ == ModuliType::P3M) ? "P3M" : "plain",
                                                            pmeOrder_,
                                                            gridSize_[XX], gridSize_[YY], gridSize_[ZZ]);
            /* This has no tester ID, but then only one tester uses testId_ currently - the correctness tester */
            std::cout << testDescription << std::endl;
            testId_ = testDescription;
        }

        //! Calls the PME init which computes the B-spline moduli
        void RunTest()
        {
            PmeInit();
        }

        //! This checks the correctness of the output
        void CheckOutputs()
        {
            gmx::test::TestReferenceData    refData;
            gmx::test::TestReferenceChecker checker(refData.rootChecker());
            checker.setDefaultTolerance(gmx::test::relativeToleranceAsUlp(1.0, 6)); // ULP == 4 failed
            for (int i = 0; i < DIM; i++)
            {
                std::string seqId = testId_ + gmx::formatString(", dimension %d", i);
                checker.checkSequenceArray(gridSize_[i], pme_.get()->bsp_mod[i], seqId.c_str());
            }
        }
};

/*! \brief Test for PME B-spline moduli computation */
TEST_P(PmeBSplineModuliTest, BehavesAsExpected)
{
    Check();
}


/* Invalid input instances */

//! Sane interpolation order
const int       sanePmeOrder = 4;
//! Sane grid size
const gmx::IVec saneGridSize = {32, 25, 47};
/*! \brief Hand-picked invalid input for the death tests */
static std::vector<BSplineModuliInputParameters> const invalidInputs
{
    /* Invalid grid sizes */
    BSplineModuliInputParameters {
        gmx::IVec {
            -1, 0, 0
        }, sanePmeOrder, ModuliType::P3M
    },
    BSplineModuliInputParameters {
        gmx::IVec {
            0, 0, 0
        }, sanePmeOrder, ModuliType::P3M
    },
    BSplineModuliInputParameters {
        gmx::IVec {
            0, 2, 0
        }, sanePmeOrder, ModuliType::PME
    },
    BSplineModuliInputParameters {
        gmx::IVec {
            INT_MAX, 0, 0
        }, sanePmeOrder, ModuliType::P3M
    },
    BSplineModuliInputParameters {
        gmx::IVec {
            64, 2, 64
        }, sanePmeOrder, ModuliType::PME
    },
    /* Invalid interpolation orders */
    BSplineModuliInputParameters {
        saneGridSize, -1, ModuliType::PME
    },
    // BSplineModuliInputParameters{saneGridSize, 0, FALSE}, //segfault
    // BSplineModuliInputParameters{saneGridSize, 1, FALSE}, //segfault
    BSplineModuliInputParameters {
        saneGridSize, INT_MAX, ModuliType::PME
    },
    BSplineModuliInputParameters {
        saneGridSize, 8 + 1, ModuliType::P3M  // P3M only supports orders up to 8
    },
    BSplineModuliInputParameters {
        saneGridSize, PME_ORDER_MAX + 1, ModuliType::PME
    },
};

/*! \brief Instantiation of the PME B-spline moduli creation test with invalid input */
INSTANTIATE_TEST_CASE_P(InsaneInput, PmeBSplineModuliTest,
                            ::testing::Combine(::testing::ValuesIn(invalidInputs),
                                                   ::testing::Values(g_TestForFatalErrorPtr)));


/* Valid input instances */

//! A couple of valid inputs - self-documenting
static std::vector<gmx::IVec> const sampleGridSizes
{
    gmx::IVec {
        64, 32, 64
    },
    gmx::IVec {
        57, 84, 29
    }
};

/*! \brief Instantiation of the PME B-spline moduli creation test with valid input - up to order of 8 */
INSTANTIATE_TEST_CASE_P(SaneInput1, PmeBSplineModuliTest, ::testing::Combine(::testing::Combine(
                                                                                         ::testing::ValuesIn(sampleGridSizes),
                                                                                         ::testing::Range(3, 8 + 1),
                                                                                         ::testing::Values(ModuliType::PME, ModuliType::P3M)
                                                                                     ), ::testing::Values(g_TestForCorrectnessPtr)));
/*! \brief Instantiation of the PME B-spline moduli creation test with valid input -
 * from order of 8 up to PME_ORDER_MAX, plain PME only */
/*
   INSTANTIATE_TEST_CASE_P(SaneInput2, PmeBSplineModuliTest, ::testing::Combine(::testing::Combine(
                                    ::testing::ValuesIn(sampleGridSizes),
                                    ::testing::Range(8 + 1, PME_ORDER_MAX + 1),
                                    ::testing::Values(ModuliType::PME)
                                ), ::testing::Values(g_TestForCorrectnessPtr)));
 */
