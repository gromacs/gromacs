/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
 * Implements PME force gathering tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include <string>

#include <gmock/gmock.h>

#include "gromacs/mdtypes/inputrec.h"
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

/* Valid input instances */

//! A couple of valid inputs for boxes.
static std::vector<Matrix3x3> const c_sampleBoxes
{
    // normal box
    Matrix3x3 {{
                   8.0f, 0.0f, 0.0f,
                   0.0f, 3.4f, 0.0f,
                   0.0f, 0.0f, 2.0f
               }},
    // triclinic box
    Matrix3x3 {{
                   7.0f, 0.0f, 0.0f,
                   0.0f, 4.1f, 0.0f,
                   3.5f, 2.0f, 12.2f
               }},
};

//! A couple of valid inputs for grid sizes
static std::vector<IVec> const c_sampleGridSizes
{
    IVec {
        16, 12, 14
    },
    IVec {
        13, 15, 11
    }
};
//! Random charges
static std::vector<real> const c_sampleChargesFull
{
    4.95f, 3.11f, 3.97f, 1.08f, 2.09f, 1.1f, 4.13f, 3.31f, 2.8f, 5.83f, 5.09f, 6.1f, 2.86f, 0.24f, 5.76f, 5.19f, 0.72f
};

//! All the input atom gridline indices
static std::vector<IVec> const c_sampleGridLineIndicesFull
{
    IVec {
        4, 2, 6
    },
    IVec {
        1, 4, 10
    },
    IVec {
        0, 6, 6
    },
    IVec {
        0, 1, 4
    },
    IVec {
        6, 3, 0
    },
    IVec {
        7, 2, 2
    },
    IVec {
        8, 3, 1
    },
    IVec {
        4, 0, 3
    },
    IVec {
        0, 0, 0
    },
    IVec {
        8, 5, 8
    },
    IVec {
        4, 4, 2
    },
    IVec {
        7, 1, 7
    },
    IVec {
        8, 5, 5
    },
    IVec {
        2, 6, 5
    },
    IVec {
        1, 6, 2
    },
    IVec {
        7, 1, 8
    },
    IVec {
        3, 5, 1
    },
};

// Spline values/derivatives below are also generated randomly, so they are bogus,
// but that should not affect the reproducibility, which we're after

//! A lot of bogus input spline values - should have at list (max PME order = 5) * (DIM = 3) * (total unique atom number in all test cases = 16) values
static std::vector<real> const c_sampleSplineValuesFull
{
    0.12f, 0.81f, 0.29f, 0.22f, 0.13f, 0.19f, 0.12f, 0.8f, 0.44f, 0.38f, 0.32f, 0.36f, 0.27f, 0.11f, 0.17f, 0.94f, 0.07f, 0.9f, 0.98f, 0.96f, 0.07f, 0.94f, 0.77f, 0.24f, 0.84f, 0.16f, 0.77f, 0.57f, 0.52f, 0.27f, 0.39f, 0.45f, 0.6f, 0.59f, 0.44f, 0.91f, 0.97f, 0.43f, 0.24f, 0.52f, 0.73f, 0.55f, 0.99f, 0.39f, 0.97f, 0.35f, 0.1f, 0.68f, 0.19f, 0.1f, 0.77f, 0.2f, 0.43f, 0.69f, 0.76f, 0.32f, 0.31f, 0.94f, 0.53f, 0.6f, 0.93f, 0.57f, 0.94f, 0.88f, 0.75f, 0.77f, 0.91f, 0.72f, 0.07f, 0.78f, 0.09f, 0.02f, 0.48f, 0.97f, 0.89f, 0.39f, 0.48f, 0.19f, 0.02f, 0.92f, 0.8f, 0.41f, 0.53f, 0.32f, 0.38f, 0.58f, 0.36f, 0.46f, 0.92f, 0.91f, 0.01f, 0.86f, 0.54f, 0.86f, 0.94f, 0.37f, 0.35f, 0.81f, 0.89f, 0.48f,
    0.34f, 0.18f, 0.11f, 0.02f, 0.87f, 0.95f, 0.66f, 0.67f, 0.38f, 0.45f, 0.04f, 0.94f, 0.54f, 0.76f, 0.58f, 0.83f, 0.31f, 0.73f, 0.71f, 0.06f, 0.35f, 0.32f, 0.35f, 0.61f, 0.27f, 0.98f, 0.83f, 0.11f, 0.3f, 0.42f, 0.95f, 0.69f, 0.58f, 0.29f, 0.1f, 0.68f, 0.94f, 0.62f, 0.51f, 0.47f, 0.04f, 0.47f, 0.34f, 0.71f, 0.52f, 0.19f, 0.69f, 0.5f, 0.59f, 0.05f, 0.74f, 0.11f, 0.4f, 0.81f, 0.24f, 0.53f, 0.71f, 0.07f, 0.17f, 0.41f, 0.23f, 0.78f, 0.27f, 0.1f, 0.71f, 0.36f, 0.67f, 0.6f, 0.94f, 0.69f, 0.19f, 0.58f, 0.68f, 0.5f, 0.62f, 0.38f, 0.29f, 0.44f, 0.04f, 0.89f, 0.0f, 0.76f, 0.22f, 0.16f, 0.08f, 0.62f, 0.51f, 0.62f, 0.83f, 0.72f, 0.96f, 0.99f, 0.4f, 0.79f, 0.83f, 0.21f, 0.43f, 0.32f, 0.44f, 0.72f,
    0.21f, 0.4f, 0.93f, 0.07f, 0.11f, 0.41f, 0.24f, 0.04f, 0.36f, 0.15f, 0.92f, 0.08f, 0.99f, 0.35f, 0.42f, 0.7f, 0.17f, 0.39f, 0.69f, 0.0f, 0.86f, 0.89f, 0.59f, 0.81f, 0.77f, 0.15f, 0.89f, 0.17f, 0.76f, 0.67f, 0.58f, 0.78f, 0.26f, 0.19f, 0.69f, 0.18f, 0.46f, 0.6f, 0.69f, 0.23f, 0.34f, 0.3f, 0.64f, 0.34f, 0.6f, 0.99f, 0.69f, 0.57f, 0.75f, 0.07f, 0.36f, 0.75f, 0.81f, 0.8f, 0.42f, 0.09f, 0.94f, 0.66f, 0.35f, 0.67f, 0.34f, 0.66f, 0.02f, 0.47f, 0.78f, 0.21f, 0.02f, 0.18f, 0.42f, 0.2f, 0.46f, 0.34f, 0.4f, 0.46f, 0.96f, 0.86f, 0.25f, 0.25f, 0.22f, 0.37f, 0.59f, 0.19f, 0.45f, 0.61f, 0.04f, 0.71f, 0.77f, 0.51f, 0.77f, 0.15f, 0.78f, 0.36f, 0.62f, 0.24f, 0.86f, 0.2f, 0.77f, 0.08f, 0.09f, 0.3f,
    0.0f, 0.6f, 0.99f, 0.69f,
};

//! A lot of bogus input spline derivatives - should have at list (max PME order = 5) * (DIM = 3) * (total unique atom number in all test cases = 16) values
static std::vector<real> const c_sampleSplineDerivativesFull
{
    0.82f, 0.88f, 0.83f, 0.11f, 0.93f, 0.32f, 0.71f, 0.37f, 0.69f, 0.88f, 0.11f, 0.38f, 0.25f, 0.5f, 0.36f, 0.81f, 0.78f, 0.31f, 0.66f, 0.32f, 0.27f, 0.35f, 0.53f, 0.83f, 0.08f, 0.08f, 0.94f, 0.71f, 0.65f, 0.24f, 0.13f, 0.01f, 0.33f, 0.65f, 0.24f, 0.53f, 0.45f, 0.84f, 0.33f, 0.97f, 0.31f, 0.7f, 0.03f, 0.31f, 0.41f, 0.76f, 0.12f, 0.3f, 0.57f, 0.65f, 0.87f, 0.99f, 0.42f, 0.97f, 0.32f, 0.39f, 0.73f, 0.23f, 0.03f, 0.67f, 0.97f, 0.57f, 0.42f, 0.38f, 0.54f, 0.17f, 0.53f, 0.54f, 0.18f, 0.8f, 0.76f, 0.13f, 0.29f, 0.83f, 0.77f, 0.56f, 0.4f, 0.87f, 0.36f, 0.18f, 0.59f, 0.04f, 0.05f, 0.61f, 0.26f, 0.91f, 0.62f, 0.16f, 0.89f, 0.23f, 0.26f, 0.59f, 0.33f, 0.2f, 0.49f, 0.41f, 0.25f, 0.4f, 0.16f, 0.83f,
    0.44f, 0.82f, 0.21f, 0.95f, 0.14f, 0.8f, 0.37f, 0.31f, 0.41f, 0.53f, 0.15f, 0.85f, 0.78f, 0.17f, 0.92f, 0.03f, 0.13f, 0.2f, 0.03f, 0.33f, 0.87f, 0.38f, 0, 0.08f, 0.79f, 0.36f, 0.53f, 0.05f, 0.07f, 0.94f, 0.23f, 0.85f, 0.13f, 0.27f, 0.23f, 0.22f, 0.26f, 0.38f, 0.15f, 0.48f, 0.18f, 0.33f, 0.23f, 0.62f, 0.1f, 0.36f, 0.99f, 0.07f, 0.02f, 0.04f, 0.09f, 0.29f, 0.52f, 0.29f, 0.83f, 0.97f, 0.61f, 0.81f, 0.49f, 0.56f, 0.08f, 0.09f, 0.03f, 0.65f, 0.46f, 0.1f, 0.06f, 0.06f, 0.39f, 0.29f, 0.04f, 0.03f, 0.1f, 0.83f, 0.94f, 0.59f, 0.97f, 0.82f, 0.2f, 0.66f, 0.23f, 0.11f, 0.03f, 0.16f, 0.27f, 0.53f, 0.94f, 0.46f, 0.43f, 0.29f, 0.97f, 0.64f, 0.46f, 0.37f, 0.43f, 0.48f, 0.37f, 0.93f, 0.5f, 0.2f,
    0.92f, 0.09f, 0.74f, 0.55f, 0.44f, 0.05f, 0.13f, 0.17f, 0.79f, 0.44f, 0.11f, 0.6f, 0.64f, 0.05f, 0.96f, 0.3f, 0.45f, 0.47f, 0.42f, 0.74f, 0.91f, 0.06f, 0.89f, 0.24f, 0.26f, 0.68f, 0.4f, 0.88f, 0.5f, 0.65f, 0.48f, 0.15f, 0.0f, 0.41f, 0.67f, 0.4f, 0.31f, 0.73f, 0.77f, 0.36f, 0.26f, 0.74f, 0.46f, 0.56f, 0.78f, 0.92f, 0.32f, 0.9f, 0.06f, 0.55f, 0.6f, 0.13f, 0.38f, 0.93f, 0.5f, 0.92f, 0.96f, 0.82f, 0.0f, 0.04f, 0.9f, 0.55f, 0.97f, 1.0f, 0.23f, 0.46f, 0.52f, 0.49f, 0.0f, 0.32f, 0.16f, 0.4f, 0.62f, 0.36f, 0.03f, 0.63f, 0.16f, 0.58f, 0.97f, 0.03f, 0.44f, 0.07f, 0.22f, 0.75f, 0.32f, 0.61f, 0.94f, 0.33f, 0.7f, 0.57f, 0.5f, 0.84f, 0.7f, 0.47f, 0.18f, 0.09f, 0.25f, 0.77f, 0.94f, 0.85f,
    0.09f, 0.83f, 0.02f, 0.91f,
};

//! 2 c_sample grids - only non-zero values have to be listed
static std::vector<SparseRealGridValuesInput> const c_sampleGrids
{
    SparseRealGridValuesInput {{
                                   IVec {
                                       0, 0, 0
                                   }, 3.5f
                               }, {
                                   IVec {
                                       7, 0, 0
                                   }, -2.5f
                               },
                               {
                                   IVec {
                                       3, 5, 7
                                   }, -0.006f
                               }, {
                                   IVec {
                                       1, 6, 7
                                   }, -2.76f
                               }, {
                                   IVec {
                                       3, 1, 2
                                   }, 0.6f
                               },  {
                                   IVec {
                                       6, 2, 4
                                   }, 7.1f
                               }, {
                                   IVec {
                                       4, 5, 6
                                   }, 4.1f
                               }, {
                                   IVec {
                                       4, 4, 6
                                   }, -3.7f
                               }, },
    SparseRealGridValuesInput {{
                                   IVec {
                                       0, 4, 0
                                   }, 6.f
                               }, {
                                   IVec {
                                       4, 2, 7
                                   }, 13.76f
                               }, {
                                   IVec {
                                       0, 6, 7
                                   }, 3.6f
                               }, {
                                   IVec {
                                       3, 2, 8
                                   }, 0.61f
                               },
                               {
                                   IVec {
                                       5, 4, 3
                                   }, 2.1f
                               },
                               {
                                   IVec {
                                       2, 5, 10
                                   }, 3.6f
                               }, {
                                   IVec {
                                       5, 3, 6
                                   }, 2.1f
                               }, {
                                   IVec {
                                       6, 6, 6
                                   }, 6.1f
                               }, }
};

//! Input forces for reduction
static std::vector<RVec> const c_sampleForcesFull {
    RVec {
        0.02f, 0.87f, 0.95f
    }, RVec {
        0.66f, 0.67f, 0.38f
    }, RVec {
        0.45f, 0.04f, 0.94f
    }, RVec {
        0.54f, 0.76f, 0.58f
    },
    RVec {
        0.83f, 0.31f, 0.73f
    }, RVec {
        0.71f, 0.06f, 0.35f
    }, RVec {
        0.32f, 0.35f, 0.61f
    }, RVec {
        0.27f, 0.98f, 0.83f
    },
    RVec {
        0.11f, 0.3f, 0.42f
    }, RVec {
        0.95f, 0.69f, 0.58f
    }, RVec {
        0.29f, 0.1f, 0.68f
    }, RVec {
        0.94f, 0.62f, 0.51f
    },
    RVec {
        0.47f, 0.04f, 0.47f
    }, RVec {
        0.34f, 0.71f, 0.52f
    }
};

//! PME orders to test
static std::vector<int> const pmeOrders {
    3, 4, 5
};
//! Atom counts to test
static std::vector<size_t> const atomCounts {
    1, 2, 13
};

/* Helper structures for test input */

//! A structure for all the spline data which depends in size both on the PME order and atom count
struct AtomAndPmeOrderSizedData
{
    //! Spline values
    SplineParamsVector splineValues;
    //! Spline derivatives
    SplineParamsVector splineDerivatives;
};

//! A structure for all the input atom data which depends in size on atom count - including a range of spline data for different PME orders
struct AtomSizedData
{
    //! Gridline indices
    GridLineIndicesVector                   gridLineIndices;
    //! Charges
    ChargesVector                           charges;
    //! Coordinates
    CoordinatesVector                       coordinates;
    //! Spline data for different orders
    std::map<int, AtomAndPmeOrderSizedData> splineDataByPmeOrder;
};

//! A range of the test input data sets, uniquely identified by the atom count
typedef std::map<size_t, AtomSizedData> InputDataByAtomCount;

/*! \brief Convenience typedef of the test input parameters - unit cell box, PME interpolation order, grid dimensions,
 *  grid values, overwriting/reducing the input forces, atom count.
 *
 * The rest of the atom-related input data - gridline indices, spline theta values, spline dtheta values, atom charges -
 * is looked up in the inputAtomDataSets_ test fixture variable.
 */
typedef std::tuple<Matrix3x3, int, IVec, SparseRealGridValuesInput, PmeForceOutputHandling, size_t> GatherInputParameters;

//! Test fixture
class PmeGatherTest : public ::testing::TestWithParam<GatherInputParameters>
{
    private:
        //! Storage of all the input atom datasets
        static InputDataByAtomCount s_inputAtomDataSets_;

    public:
        //! Default constructor
        PmeGatherTest()  = default;
        //! Sets the input atom data references once
        static void SetUpTestCase()
        {
            auto                gridLineIndicesIt = c_sampleGridLineIndicesFull.begin();
            auto                chargesIt         = c_sampleChargesFull.begin();
            for (auto atomCount : atomCounts)
            {
                AtomSizedData atomData;
                atomData.gridLineIndices = GridLineIndicesVector::fromVector(gridLineIndicesIt, gridLineIndicesIt + atomCount);
                gridLineIndicesIt       += atomCount;
                atomData.charges         = ChargesVector::fromVector(chargesIt, chargesIt + atomCount);
                chargesIt               += atomCount;
                atomData.coordinates.resize(atomCount, RVec {1e6, 1e7, -1e8});
                /* The coordinates are intentionally bogus in this test - only the size matters; the gridline indices are fed directly as inputs */
                for (auto pmeOrder : pmeOrders)
                {
                    AtomAndPmeOrderSizedData splineData;
                    const size_t             dimSize = atomCount * pmeOrder;
                    for (int dimIndex = 0; dimIndex < DIM; dimIndex++)
                    {
                        splineData.splineValues[dimIndex] = SplineParamsDimVector::fromVector(c_sampleSplineValuesFull.begin() + dimIndex * dimSize,
                                                                                              c_sampleSplineValuesFull.begin() + (dimIndex + 1) * dimSize);
                        splineData.splineDerivatives[dimIndex] = SplineParamsDimVector::fromVector(c_sampleSplineDerivativesFull.begin() + dimIndex * dimSize,
                                                                                                   c_sampleSplineDerivativesFull.begin() + (dimIndex + 1) * dimSize);
                    }
                    atomData.splineDataByPmeOrder[pmeOrder] = splineData;
                }
                s_inputAtomDataSets_[atomCount] = atomData;
            }
        }
        //! The test
        void runTest()
        {
            /* Getting the input */
            Matrix3x3                 box;
            int                       pmeOrder;
            IVec                      gridSize;
            size_t                    atomCount;
            SparseRealGridValuesInput nonZeroGridValues;
            PmeForceOutputHandling    inputForceTreatment;
            std::tie(box, pmeOrder, gridSize, nonZeroGridValues, inputForceTreatment, atomCount) = GetParam();
            auto inputAtomData       = s_inputAtomDataSets_[atomCount];
            auto inputAtomSplineData = inputAtomData.splineDataByPmeOrder[pmeOrder];

            /* Storing the input where it's needed, running the test */
            t_inputrec inputRec;
            inputRec.nkx         = gridSize[XX];
            inputRec.nky         = gridSize[YY];
            inputRec.nkz         = gridSize[ZZ];
            inputRec.pme_order   = pmeOrder;
            inputRec.coulombtype = eelPME;
            inputRec.epsilon_r   = 1.0;

            TestReferenceData                     refData;
            const std::map<CodePath, std::string> modesToTest = {{CodePath::CPU, "CPU"},
                                                                 {CodePath::CUDA, "CUDA"}};
            for (const auto &mode : modesToTest)
            {
                const bool supportedInput = pmeSupportsInputForMode(&inputRec, mode.first);
                if (!supportedInput)
                {
                    /* Testing the failure for the unsupported input */
                    EXPECT_THROW(pmeInitAtoms(&inputRec, mode.first, nullptr, inputAtomData.coordinates, inputAtomData.charges, box), NotImplementedError);
                    continue;
                }

                const auto contextsToTest = pmeEnv->getHardwareContexts(mode.first);
                for (const auto &context : contextsToTest)
                {
                    /* Describing the test uniquely */
                    SCOPED_TRACE(formatString("Testing force gathering with %s %sfor PME grid size %d %d %d"
                                              ", order %d, %zu atoms, %s",
                                              mode.second.c_str(), context.getDescription().c_str(),
                                              gridSize[XX], gridSize[YY], gridSize[ZZ],
                                              pmeOrder,
                                              atomCount,
                                              (inputForceTreatment == PmeForceOutputHandling::ReduceWithInput) ? "with reduction" : "without reduction"
                                              ));

                    PmeSafePointer pmeSafe = pmeInitAtoms(&inputRec, mode.first, context.getDeviceInfo(), inputAtomData.coordinates, inputAtomData.charges, box);

                    /* Setting some more inputs */
                    pmeSetRealGrid(pmeSafe.get(), mode.first, nonZeroGridValues);
                    pmeSetGridLineIndices(pmeSafe.get(), mode.first, inputAtomData.gridLineIndices);
                    for (int dimIndex = 0; dimIndex < DIM; dimIndex++)
                    {
                        pmeSetSplineData(pmeSafe.get(), mode.first, inputAtomSplineData.splineValues[dimIndex], PmeSplineDataType::Values, dimIndex);
                        pmeSetSplineData(pmeSafe.get(), mode.first, inputAtomSplineData.splineDerivatives[dimIndex], PmeSplineDataType::Derivatives, dimIndex);
                    }

                    /* Explicitly copying the sample forces to be able to modify them */
                    auto inputForcesFull(c_sampleForcesFull);
                    GMX_RELEASE_ASSERT(inputForcesFull.size() >= atomCount, "Bad input forces size");
                    auto forces = ForcesVector::fromVector(inputForcesFull.begin(), inputForcesFull.begin() + atomCount);

                    /* Running the force gathering itself */
                    pmePerformGather(pmeSafe.get(), mode.first, inputForceTreatment, forces);
                    pmeFinalizeTest(pmeSafe.get(), mode.first);

                    /* Check the output forces correctness */
                    TestReferenceChecker forceChecker(refData.rootChecker());
                    const auto           ulpTolerance = 3 * pmeOrder;
                    forceChecker.setDefaultTolerance(relativeToleranceAsUlp(1.0, ulpTolerance));
                    forceChecker.checkSequence(forces.begin(), forces.end(), "Forces");
                }
            }
        }
};

// An instance of static atom data
InputDataByAtomCount PmeGatherTest::s_inputAtomDataSets_;

//! Test for PME force gathering
TEST_P(PmeGatherTest, ReproducesOutputs)
{
    EXPECT_NO_THROW(runTest());
}

//! Instantiation of the PME gathering test
INSTANTIATE_TEST_CASE_P(SaneInput, PmeGatherTest, ::testing::Combine(::testing::ValuesIn(c_sampleBoxes),
                                                                         ::testing::ValuesIn(pmeOrders),
                                                                         ::testing::ValuesIn(c_sampleGridSizes),
                                                                         ::testing::ValuesIn(c_sampleGrids),
                                                                         ::testing::Values(PmeForceOutputHandling::Set, PmeForceOutputHandling::ReduceWithInput),
                                                                         ::testing::ValuesIn(atomCounts)));

}
}
}
