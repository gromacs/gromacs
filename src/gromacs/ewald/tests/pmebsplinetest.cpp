/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Implements PME B-spline moduli computation tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include <array>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_internal.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "pmetestcommon.h"

namespace gmx
{
namespace test
{
namespace
{
/*! \brief Moduli algorithm type */
enum class ModuliType
{
    PME,
    P3M
};

/*! \brief Convenience typedef of input parameters - grid dimensions, PME interpolation order, moduli type */
typedef std::tuple<IVec, int, ModuliType> BSplineModuliInputParameters;

/*! \brief Test fixture for testing PME B-spline moduli creation */
class PmeBSplineModuliTest : public ::testing::TestWithParam<BSplineModuliInputParameters>
{
public:
    PmeBSplineModuliTest() = default;
    //! The whole logic being tested is contained here
    static void runTest()
    {
        /* Getting the input */
        const BSplineModuliInputParameters parameters = GetParam();
        int                                pmeOrder;
        IVec                               gridSize;
        ModuliType                         moduliType;
        std::tie(gridSize, pmeOrder, moduliType) = parameters;

        /* Describing the test in case it fails */
        SCOPED_TRACE(formatString(
                "Testing B-spline moduli creation (%s) for PME order %d, grid size %d %d %d",
                (moduliType == ModuliType::P3M) ? "P3M" : "plain",
                pmeOrder,
                gridSize[XX],
                gridSize[YY],
                gridSize[ZZ]));

        /* Storing the input where it's needed */
        t_inputrec inputRec;
        inputRec.nkx         = gridSize[XX];
        inputRec.nky         = gridSize[YY];
        inputRec.nkz         = gridSize[ZZ];
        inputRec.coulombtype = (moduliType == ModuliType::P3M) ? CoulombInteractionType::P3mAD
                                                               : CoulombInteractionType::Pme;
        inputRec.pme_order   = pmeOrder;

        /* PME initialization call which checks the inputs and computes the B-spline moduli according to the grid sizes. */
        PmeSafePointer pme = pmeInitEmpty(&inputRec);

        /* Setting up the checker */
        TestReferenceData    refData;
        TestReferenceChecker checker(refData.rootChecker());
        checker.setDefaultTolerance(relativeToleranceAsPrecisionDependentUlp(
                1.0, c_splineModuliSinglePrecisionUlps, getSplineModuliDoublePrecisionUlps(pmeOrder + 1)));

        /* Perform a correctness check */
        const char* dimString[] = { "X", "Y", "Z" };
        for (int i = 0; i < DIM; i++)
        {
            checker.checkSequenceArray(gridSize[i], pme->bsp_mod[i].data(), dimString[i]);
        }
    }
};

// Different aliases are needed to reuse the test class without instantiating tests more than once
//! Failure test alias
using PmeBSplineModuliFailureTest = PmeBSplineModuliTest;
TEST_P(PmeBSplineModuliFailureTest, Throws)
{
    EXPECT_THROW_GMX(runTest(), InconsistentInputError);
}
//! Correctness test alias
using PmeBSplineModuliCorrectnessTest = PmeBSplineModuliTest;
TEST_P(PmeBSplineModuliCorrectnessTest, ReproducesValues)
{
    EXPECT_NO_THROW_GMX(runTest());
}

/* Invalid input instances */

//! Sane interpolation order
const int sanePmeOrder = 4;
//! Sane grid size
const IVec saneGridSize = { 32, 25, 47 };
/*! \brief Hand-picked invalid input for the exception tests */
std::vector<BSplineModuliInputParameters> const invalidInputs{
    /* Invalid grid sizes */
    BSplineModuliInputParameters{ IVec{ -1, 10, 10 }, sanePmeOrder, ModuliType::P3M },
    BSplineModuliInputParameters{ IVec{ 40, 0, 20 }, sanePmeOrder, ModuliType::P3M },
    BSplineModuliInputParameters{ IVec{ 64, 2, 64 }, sanePmeOrder, ModuliType::PME },
    /* Invalid interpolation orders */
    BSplineModuliInputParameters{
            saneGridSize,
            8 + 1,
            ModuliType::P3M // P3M only supports orders up to 8
    },
    BSplineModuliInputParameters{ saneGridSize, PME_ORDER_MAX + 1, ModuliType::PME },
};

/*! \brief Instantiation of the PME B-spline moduli creation test with invalid input */
INSTANTIATE_TEST_SUITE_P(InsaneInput, PmeBSplineModuliFailureTest, ::testing::ValuesIn(invalidInputs));

/* Valid input instances */

//! A couple of valid inputs for grid sizes. It is good to test both even and odd dimensions.
std::vector<IVec> const sampleGridSizes{ IVec{ 64, 32, 64 }, IVec{ 57, 84, 29 } };

/*! \brief Instantiation of the PME B-spline moduli creation test with valid input - up to order of 8 */
INSTANTIATE_TEST_SUITE_P(SaneInput1,
                         PmeBSplineModuliCorrectnessTest,
                         ::testing::Combine(::testing::ValuesIn(sampleGridSizes),
                                            ::testing::Range(3, 8 + 1),
                                            ::testing::Values(ModuliType::PME, ModuliType::P3M)));
} // namespace
} // namespace test
} // namespace gmx
