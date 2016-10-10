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

#include <string>

#include <gmock/gmock.h>

#include "gromacs/ewald/pme-internal.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "pmebasetester.h"

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
class PmeBSplineModuliTest : public ::testing::TestWithParam<BSplineModuliTestParameters>
{
    public:
        PmeBSplineModuliTest() {}
};

//! The whole logic being tested is contained here
void RunTest(const BSplineModuliInputParameters &parameters)
{
    /* Getting the input */
    int         pmeOrder;
    gmx::IVec   gridSize;
    ModuliType  moduliType;
    std::tie(gridSize, pmeOrder, moduliType) = parameters;

    /* Storing the input where it's needed */
    t_inputrec *inputRec = s_mdModules.inputrec();
    inputRec->nkx         = gridSize[XX];
    inputRec->nky         = gridSize[YY];
    inputRec->nkz         = gridSize[ZZ];
    inputRec->coulombtype = (moduliType == ModuliType::P3M) ? eelP3M_AD : eelPME;
    inputRec->pme_order   = pmeOrder;

    /* Describing the test uniquely */
    std::string testDescription = gmx::formatString("Testing B-spline moduli creation (%s) for PME order %d, grid size %d %d %d",
                                                    (moduliType == ModuliType::P3M) ? "P3M" : "plain",
                                                    pmeOrder,
                                                    gridSize[XX], gridSize[YY], gridSize[ZZ]);
    std::cout << testDescription << std::endl;
    std::string testId = testDescription;

    /* This computes the modules */
    pmeSafePointer pme = PmeDummyInit(inputRec);

    /* Do a correctness check unless we are dead after the PME init */
    gmx::test::TestReferenceData    refData;
    gmx::test::TestReferenceChecker checker(refData.rootChecker());

    gmx_uint64_t                    ulpTolerance = (GMX_DOUBLE && (moduliType == ModuliType::P3M)) ? 10 : 6;
    /* P3M moduli use polynomial expansions that are strongly affected by rounding errors. */
    checker.setDefaultTolerance(gmx::test::relativeToleranceAsUlp(1.0, ulpTolerance));
    for (int i = 0; i < DIM; i++)
    {
        std::string seqId = testId + gmx::formatString(", dimension %d", i);
        checker.checkSequenceArray(gridSize[i], pme.get()->bsp_mod[i], seqId.c_str());
    }
}

/*! \brief Test for PME B-spline moduli computation */
TEST_P(PmeBSplineModuliTest, BehavesAsExpected)
{
    BSplineModuliInputParameters parameters;
    TesterFuncPtr                testerFunction;
    std::tie(parameters, testerFunction) = GetParam();
    GMX_ASSERT((void *)testerFunction != NULL, "No tester function specified");

    TestFunc test = std::bind(RunTest, parameters);
    testerFunction(test);
}


/* Invalid input instances */

//! Sane interpolation order
const int       sanePmeOrder = 4;
//! Sane grid size
const gmx::IVec saneGridSize = {32, 25, 47};
/*! \brief Hand-picked invalid input for the exception tests */
static std::vector<BSplineModuliInputParameters> const invalidInputs
{
    /* Invalid grid sizes */
    BSplineModuliInputParameters {
        gmx::IVec {
            -1, 10, 10
        }, sanePmeOrder, ModuliType::P3M
    },
    BSplineModuliInputParameters {
        gmx::IVec {
            40, 0, 20
        }, sanePmeOrder, ModuliType::P3M
    },
    BSplineModuliInputParameters {
        gmx::IVec {
            64, 2, 64
        }, sanePmeOrder, ModuliType::PME
    },
    /* Invalid interpolation orders */
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
                                                   ::testing::Values(g_TestForInvalidInputExceptionPtr)));

/* Valid input instances */

//! A couple of valid inputs for grid sizes. It is good to test both even and odd dimensions.
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
                                                                                     ), ::testing::Values(g_TestForNoExceptionPtr)));
/*! \brief Instantiation of the PME B-spline moduli creation test with valid input -
 * from order of 8 up to PME_ORDER_MAX, plain PME only */
/*
   INSTANTIATE_TEST_CASE_P(SaneInput2, PmeBSplineModuliTest, ::testing::Combine(::testing::Combine(
                                    ::testing::ValuesIn(sampleGridSizes),
                                    ::testing::Range(8 + 1, PME_ORDER_MAX + 1),
                                    ::testing::Values(ModuliType::PME)
                                ), ::testing::Values(g_TestForNoExceptionPtr)));
 */
