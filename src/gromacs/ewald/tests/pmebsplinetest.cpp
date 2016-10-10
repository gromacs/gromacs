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
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/scoped_cptr.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"

#include "pmeabstracttest.h"

/*! \brief Moduli algorithm type */
enum class ModuliType
{
    PME,
    P3M
};
// When adding new members to enums, don't forget to alter the test instantiations!

/*! \brief Convenience typedef of test parameters - grid dimensions (ivec), PME interpolation order, moduli type */
typedef std::tuple<const int *, int, ModuliType> BSplineModuliTestParameters;

/*! \brief Abstract test fixture for testing PME B-spline moduli creation */
class AbstractPmeBSplineModuliTestBase : virtual public AbstractPmeTestBase,
                                         virtual public ::testing::TestWithParam<BSplineModuliTestParameters>
{
    public:
        //! PME spline interpolation order
        int         pmeOrder_;
        //! PME grid size - is actualy an ivec
        const int  *gridSize_;
        //! Moduli algorithm type - see above
        ModuliType  moduliType_;
        //! Unique test ID for the reference data
        std::string testId_;

        /*! \brief An initializer */
        void SetUp()
        {
            /* Getting the input */
            std::tie(gridSize_, pmeOrder_, moduliType_) = GetParam();

            /* Storing the input where it's needed */
            inputRec_->nkx         = gridSize_[XX];
            inputRec_->nky         = gridSize_[YY];
            inputRec_->nkz         = gridSize_[ZZ];
            inputRec_->coulombtype = (moduliType_ == ModuliType::P3M) ? eelP3M_AD : eelPME;
            inputRec_->pme_order   = pmeOrder_;

            /* Describing the test for output */
            std::string testDescription = gmx::formatString("Testing B-spline moduli creation (%s) for PME order %d, grid size %d %d %d",
                                                            (moduliType_ == ModuliType::P3M) ? "P3M" : "plain",
                                                            pmeOrder_,
                                                            gridSize_[XX], gridSize_[YY], gridSize_[ZZ]);
            std::cout << testDescription << std::endl;
            testId_ = testDescription;
        }
};


/* Invalid input death test */

/*! \brief Test fixture for testing PME B-spline moduli creation termination on invalid input */
class PmeBSplineModuliDeathTest : virtual public AbstractPmeBSplineModuliTestBase, public AbstractPmeDeathTestBase
{
};

/*! \brief Test or testing PME B-spline moduli creation termination on invalid input */
TEST_P(PmeBSplineModuliDeathTest, Terminates)
{
    Check();
}

//! Self-documenting
const int  sanePmeOrder = 4;
//! Self-documenting
const ivec saneGridSize = {32, 25, 47};
/*! \brief Hand-picked invalid input for the death tests */
static std::vector<BSplineModuliTestParameters> const invalidInputs
{
    /* Invalid grid sizes */
    BSplineModuliTestParameters {
        ivec {
            -1, 0, 0
        }, sanePmeOrder, ModuliType::P3M
    },
    BSplineModuliTestParameters {
        ivec {
            0, 0, 0
        }, sanePmeOrder, ModuliType::P3M
    },
    BSplineModuliTestParameters {
        ivec {
            0, 2, 0
        }, sanePmeOrder, ModuliType::PME
    },
    BSplineModuliTestParameters {
        ivec {
            INT_MAX, 0, 0
        }, sanePmeOrder, ModuliType::P3M
    },
    BSplineModuliTestParameters {
        ivec {
            64, 2, 64
        }, sanePmeOrder, ModuliType::PME
    },
    /* Invalid interpolation orders */
    BSplineModuliTestParameters {
        saneGridSize, -1, ModuliType::PME
    },
    // BSplineModuliTestParameters{saneGridSize, 0, FALSE}, //segfault
    // BSplineModuliTestParameters{saneGridSize, 1, FALSE}, //segfault
    BSplineModuliTestParameters {
        saneGridSize, INT_MAX, ModuliType::PME
    },
    BSplineModuliTestParameters {
        saneGridSize, 8 + 1, ModuliType::P3M  // P3M only supports orders up to 8
    },
    BSplineModuliTestParameters {
        saneGridSize, PME_ORDER_MAX + 1, ModuliType::PME
    },
};

/*! \brief Instantiation of the PME B-spline moduli creation test with invalid input */
INSTANTIATE_TEST_CASE_P(InsaneInput, PmeBSplineModuliDeathTest, ::testing::ValuesIn(invalidInputs));


/* Valid input test */

/*! \brief Test fixture for checking PME B-spline moduli creation output on valid input */
class PmeBSplineModuliCorrectnessTest : virtual public AbstractPmeBSplineModuliTestBase, public AbstractPmeCorrectnessTestBase
{
    //! This checks the correctness of the output
    void CheckOutputs()
    {
        //TODO checker_.setDefaultTolerance(
        //       gmx::test::relativeToleranceAsPrecisionDependentUlp(10.0, 64, 512));
        for (int i = 0; i < DIM; i++)
        {
            std::string seqId = testId_ + gmx::formatString(", dimension %d", i);
            checker_.checkSequenceArray(gridSize_[i], pme_.get()->bsp_mod[i], seqId.c_str());
        }
    }
};

/*! \brief Test for checking PME B-spline moduli creation on valid input */
TEST_P(PmeBSplineModuliCorrectnessTest, ReproducesBSplineModuli)
{
    Check();
}

//! A couple of valid inputs - self-documenting
static std::vector<const int *> const sampleGridSizes
{
    ivec {
        64, 32, 64
    },
    ivec {
        57, 84, 29
    }
};

/*! \brief Instantiation of the PME B-spline moduli creation test with valid input - up to order of 8 */
INSTANTIATE_TEST_CASE_P(SaneInput1, PmeBSplineModuliCorrectnessTest, ::testing::Combine(
                                    ::testing::ValuesIn(sampleGridSizes),
                                    ::testing::Range(3, 8 + 1),
                                    ::testing::Values(ModuliType::PME, ModuliType::P3M)
                                ));
/*! \brief Instantiation of the PME B-spline moduli creation test with valid input - from order of 8 up to PME_ORDER_MAX, plain PME only */
INSTANTIATE_TEST_CASE_P(SaneInput2, PmeBSplineModuliCorrectnessTest, ::testing::Combine(
                                    ::testing::ValuesIn(sampleGridSizes),
                                    ::testing::Range(8 + 1, PME_ORDER_MAX + 1),
                                    ::testing::Values(ModuliType::PME)
                                ));
