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
 * Implements PME spline computation and charge spreading tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include <string>

#include <gmock/gmock.h>

#include "gromacs/ewald/pme-internal.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "pmetestcommon.h"

/*! \brief Convenience typedef of input parameters - unit cell box, PME interpolation order, grid dimensions,
 * particle coordinates, particle charges, PME code path being tested
 * TODO: consider inclusion of local grid offsets/sizes or PME nodes counts to test the PME DD
 */
typedef std::tuple<gmx::Matrix3x3, int, gmx::IVec, CoordinatesVector, ChargesVector, PmeCodePath> SplineSpreadInputParameters;

//! Test fixture
class PmeSplineSpreadTest : public ::testing::TestWithParam<SplineSpreadInputParameters>
{
    public:
        //! Default constructor
        PmeSplineSpreadTest() = default;
        //! The test
        void RunTest()
        {
            /* Getting the input */
            gmx::Matrix3x3          box;
            int                     pmeOrder;
            gmx::IVec               gridSize;
            CoordinatesVector       coordinates;
            ChargesVector           charges;
            PmeCodePath             mode;

            std::tie(box, pmeOrder, gridSize, coordinates, charges, mode) = GetParam();
            const size_t atomCount = coordinates.size();

            /* Describing the test uniquely in case it fails */
            SCOPED_TRACE(gmx::formatString("Testing spline computation and spreading for PME grid size %d %d %d"
                                           ", order %d, %zu atoms",
                                           gridSize[XX], gridSize[YY], gridSize[ZZ],
                                           pmeOrder,
                                           atomCount));

            /* Storing the input where it's needed */
            t_inputrec *inputRec  = s_mdModules.inputrec();
            inputRec->nkx         = gridSize[XX];
            inputRec->nky         = gridSize[YY];
            inputRec->nkz         = gridSize[ZZ];
            inputRec->pme_order   = pmeOrder;
            inputRec->coulombtype = eelPME;

            pmeSafePointer  pmeSafe = PmeInitWithAtoms(inputRec, coordinates, charges, box);
            gmx_pme_t      *pme     = pmeSafe.get();
            pme_atomcomm_t *atc     = &(pmeSafe.get()->atc[0]);

            /* Running the test */
            const bool computeSplines = true, spreadCharges = true; //move?
            PmePerformSpread(pmeSafe, mode, computeSplines, spreadCharges);

            /* Outputs correctness check */
            const char                     *dimString[] = { "X", "Y", "Z" };
            gmx::test::TestReferenceData    refData;
            gmx::test::TestReferenceChecker checker(refData.rootChecker());

            if (computeSplines)
            {
                gmx_uint64_t                    ulpTolerance = 24; // empiric - quite large, no? fix for doubles!
                checker.setDefaultTolerance(gmx::test::relativeToleranceAsUlp(1.0, ulpTolerance));

                const size_t                    threadIndex = 0; // relying on running single threaded on CPU
                for (int i = 0; i < DIM; i++)
                {
                    /* Particle spline values */
                    std::string seqId = gmx::formatString("Spline values, %s", dimString[i]);
                    checker.checkSequenceArray(atomCount * pmeOrder, atc->spline[threadIndex].theta[i], seqId.c_str());
                    /* Particle spline derivatives */
                    seqId = gmx::formatString("Spline derivatives, %s", dimString[i]);
                    checker.checkSequenceArray(atomCount * pmeOrder, atc->spline[threadIndex].dtheta[i], seqId.c_str());
                }
                /* Fractional particle coordinates */
                checker.checkSequenceArray(atomCount, atc->fractx, "Fractional coordinates");
                /* Particle gridline indices */
                checker.checkSequenceArray(atomCount, atc->idx, "Gridline indices");
            }

            if (spreadCharges)
            {
                const size_t gridIndex = 0;
                gmx::IVec    localFftNdataUnused, localFftOffsetUnused, localFftSize;
                gmx_parallel_3dfft_real_limits(pme->pfft_setup[gridIndex],
                                               localFftNdataUnused,
                                               localFftOffsetUnused,
                                               localFftSize);

                gmx_uint64_t                    ulpTolerance = 72; //!
                checker.setDefaultTolerance(gmx::test::relativeToleranceAsUlp(1.0, ulpTolerance));

                /* Spread's only output is the wrapped grid */

                gmx::test::TestReferenceChecker gridValuesChecker(checker.checkCompound("NonZeroGridValues", "RealGrid"));

                for (int ix = 0; ix < gridSize[XX]; ix++)
                {
                    for (int iy = 0; iy < gridSize[YY]; iy++)
                    {
                        for (int iz = 0; iz < gridSize[ZZ]; iz++)
                        {
                            std::string  valueId        = gmx::formatString("Cell %d %d %d", ix, iy, iz);
                            const size_t gridValueIndex = (ix * localFftSize[YY] + iy) * localFftSize[ZZ] + iz;
                            const real   value          = pme->fftgrid[gridIndex][gridValueIndex];
                            const bool   empty          = (value == 0.0);
                            gridValuesChecker.checkPresent(!empty, valueId.c_str());
                            if (!empty)
                            {
                                gridValuesChecker.checkReal(value, valueId.c_str());
                            }
                        }
                    }
                }
            }
        }
};


/*! \brief Test for PME B-spline moduli computation */
TEST_P(PmeSplineSpreadTest, ReproducesOutputs)
{
    EXPECT_NO_THROW(RunTest());
}

/* Valid input instances */
//! A couple of valid inputs for boxes.
static std::vector<gmx::Matrix3x3> const sampleBoxes
{
    // normal box
    gmx::Matrix3x3 {
        8.0f, 0.0f, 0.0f,
        0.0f, 3.4f, 0.0f,
        0.0f, 0.0f, 2.0f
    },
    // triclinic box
    gmx::Matrix3x3 {
        7.0f, 0.0f, 0.0f,
        0.0f, 2.1f, 0.0f,
        6.0f, 5.0f, 3.2f
    },
};

//! A couple of valid inputs for grid sizes.
static std::vector<gmx::IVec> const sampleGridSizes
{
    gmx::IVec {
        16, 12, 14
    },
    gmx::IVec {
        9, 7, 11
    }
};

//! Random charges
static ChargesVector const sampleChargesFull
{
    4.95f, 3.11f, 3.97f, 1.08f, 2.09f, 1.1f, 4.13f, 3.31f, 2.8f, 5.83f, 5.09f, 6.1f, 2.86f, 0.24f, 5.76f, 5.19f, 0.72f
};
//! 1 charge
static ChargesVector const sampleCharges1(sampleChargesFull.begin(), sampleChargesFull.begin() + 1);
//! 2 charges
static ChargesVector const sampleCharges2(sampleChargesFull.begin() + 1, sampleChargesFull.begin() + 3);
//! 13 charges
static ChargesVector const sampleCharges13(sampleChargesFull.begin() + 3, sampleChargesFull.begin() + 16);

//! Random coordinate vectors
static CoordinatesVector const sampleCoordinatesFull
{
    {
        1.59f, 1.37f, 0.95f
    }, {
        0.03f, 1.6f, 0.22f
    }, {
        0.33f, 0.92f, 1.56f
    }, {
        1.16f, 0.75f, 0.39f
    }, {
        0.79f, 1.6f, 1.12f
    }, {
        0.49f, 1.02f, 0.22f
    }, {
        0.65f, 1.52f, 1.19f
    }, {
        1.43f, 1.1f, 0.08f
    }, {
        1.08f, 1.19f, 0.08f
    }, {
        1.6f, 0.93f, 0.53f
    }, {
        1.32f, 1.48f, 0.16f
    }, {
        0.87f, 0.22f, 0.33f
    }, {
        0.95f, 0.13f, 0.48f
    }, {
        1.23f, 0.91f, 0.68f
    }, {
        0.19f, 1.45f, 0.94f
    }, {
        1.28f, 0.46f, 0.38f
    }, {
        1.21f, 0.23f, 1.0f
    }
};
//! 1 coordinate vector
static CoordinatesVector const sampleCoordinates1(sampleCoordinatesFull.begin(), sampleCoordinatesFull.begin() + 1);
//! 2 coordinate vectors
static CoordinatesVector const sampleCoordinates2(sampleCoordinatesFull.begin() + 1, sampleCoordinatesFull.begin() + 3);
//! 13 coordinate vectors
static CoordinatesVector const sampleCoordinates13(sampleCoordinatesFull.begin() + 3, sampleCoordinatesFull.begin() + 16);

/*! \brief Instantiation of the PME spline computation test with valid input and 1 atom */
INSTANTIATE_TEST_CASE_P(SaneInput1, PmeSplineSpreadTest, ::testing::Combine(
                                    ::testing::ValuesIn(sampleBoxes),
                                    ::testing::Range(4, 5 + 1),
                                    ::testing::ValuesIn(sampleGridSizes),
                                    ::testing::Values(sampleCoordinates1),
                                    ::testing::Values(sampleCharges1),
                                    ::testing::Values(PmeCodePath::CPU)
                                ));
/*! \brief Instantiation of the PME spline computation test with valid input and 2 atoms */
INSTANTIATE_TEST_CASE_P(SaneInput2, PmeSplineSpreadTest, ::testing::Combine(
                                    ::testing::ValuesIn(sampleBoxes),
                                    ::testing::Range(4, 5 + 1),
                                    ::testing::ValuesIn(sampleGridSizes),
                                    ::testing::Values(sampleCoordinates2),
                                    ::testing::Values(sampleCharges2),
                                    ::testing::Values(PmeCodePath::CPU)
                                ));
/*! \brief Instantiation of the PME spline computation test with valid input and 13 atoms */
INSTANTIATE_TEST_CASE_P(SaneInput13, PmeSplineSpreadTest, ::testing::Combine(
                                    ::testing::ValuesIn(sampleBoxes),
                                    ::testing::Range(4, 5 + 1),
                                    ::testing::ValuesIn(sampleGridSizes),
                                    ::testing::Values(sampleCoordinates13),
                                    ::testing::Values(sampleCharges13),
                                    ::testing::Values(PmeCodePath::CPU)
                                ));
