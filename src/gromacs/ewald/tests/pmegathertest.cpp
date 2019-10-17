/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
std::vector<Matrix3x3> const c_sampleBoxes{
    // normal box
    Matrix3x3{ { 8.0F, 0.0F, 0.0F, 0.0F, 3.4F, 0.0F, 0.0F, 0.0F, 2.0F } },
    // triclinic box
    Matrix3x3{ { 7.0F, 0.0F, 0.0F, 0.0F, 4.1F, 0.0F, 3.5F, 2.0F, 12.2F } },
};

//! A couple of valid inputs for grid sizes
std::vector<IVec> const c_sampleGridSizes{ IVec{ 16, 12, 14 }, IVec{ 13, 15, 11 } };
//! Random charges
std::vector<real> const c_sampleChargesFull{ 4.95F, 3.11F, 3.97F, 1.08F, 2.09F, 1.1F,
                                             4.13F, 3.31F, 2.8F,  5.83F, 5.09F, 6.1F,
                                             2.86F, 0.24F, 5.76F, 5.19F, 0.72F };

//! All the input atom gridline indices
std::vector<IVec> const c_sampleGridLineIndicesFull{
    IVec{ 4, 2, 6 }, IVec{ 1, 4, 10 }, IVec{ 0, 6, 6 }, IVec{ 0, 1, 4 }, IVec{ 6, 3, 0 },
    IVec{ 7, 2, 2 }, IVec{ 8, 3, 1 },  IVec{ 4, 0, 3 }, IVec{ 0, 0, 0 }, IVec{ 8, 5, 8 },
    IVec{ 4, 4, 2 }, IVec{ 7, 1, 7 },  IVec{ 8, 5, 5 }, IVec{ 2, 6, 5 }, IVec{ 1, 6, 2 },
    IVec{ 7, 1, 8 }, IVec{ 3, 5, 1 },
};

// Spline values/derivatives below are also generated randomly, so they are bogus,
// but that should not affect the reproducibility, which we're after

//! A lot of bogus input spline values - should have at list (max PME order = 5) * (DIM = 3) * (total unique atom number in all test cases = 16) values
std::vector<real> const c_sampleSplineValuesFull{
    0.12F, 0.81F, 0.29F, 0.22F, 0.13F, 0.19F, 0.12F, 0.8F,  0.44F, 0.38F, 0.32F, 0.36F, 0.27F,
    0.11F, 0.17F, 0.94F, 0.07F, 0.9F,  0.98F, 0.96F, 0.07F, 0.94F, 0.77F, 0.24F, 0.84F, 0.16F,
    0.77F, 0.57F, 0.52F, 0.27F, 0.39F, 0.45F, 0.6F,  0.59F, 0.44F, 0.91F, 0.97F, 0.43F, 0.24F,
    0.52F, 0.73F, 0.55F, 0.99F, 0.39F, 0.97F, 0.35F, 0.1F,  0.68F, 0.19F, 0.1F,  0.77F, 0.2F,
    0.43F, 0.69F, 0.76F, 0.32F, 0.31F, 0.94F, 0.53F, 0.6F,  0.93F, 0.57F, 0.94F, 0.88F, 0.75F,
    0.77F, 0.91F, 0.72F, 0.07F, 0.78F, 0.09F, 0.02F, 0.48F, 0.97F, 0.89F, 0.39F, 0.48F, 0.19F,
    0.02F, 0.92F, 0.8F,  0.41F, 0.53F, 0.32F, 0.38F, 0.58F, 0.36F, 0.46F, 0.92F, 0.91F, 0.01F,
    0.86F, 0.54F, 0.86F, 0.94F, 0.37F, 0.35F, 0.81F, 0.89F, 0.48F, 0.34F, 0.18F, 0.11F, 0.02F,
    0.87F, 0.95F, 0.66F, 0.67F, 0.38F, 0.45F, 0.04F, 0.94F, 0.54F, 0.76F, 0.58F, 0.83F, 0.31F,
    0.73F, 0.71F, 0.06F, 0.35F, 0.32F, 0.35F, 0.61F, 0.27F, 0.98F, 0.83F, 0.11F, 0.3F,  0.42F,
    0.95F, 0.69F, 0.58F, 0.29F, 0.1F,  0.68F, 0.94F, 0.62F, 0.51F, 0.47F, 0.04F, 0.47F, 0.34F,
    0.71F, 0.52F, 0.19F, 0.69F, 0.5F,  0.59F, 0.05F, 0.74F, 0.11F, 0.4F,  0.81F, 0.24F, 0.53F,
    0.71F, 0.07F, 0.17F, 0.41F, 0.23F, 0.78F, 0.27F, 0.1F,  0.71F, 0.36F, 0.67F, 0.6F,  0.94F,
    0.69F, 0.19F, 0.58F, 0.68F, 0.5F,  0.62F, 0.38F, 0.29F, 0.44F, 0.04F, 0.89F, 0.0F,  0.76F,
    0.22F, 0.16F, 0.08F, 0.62F, 0.51F, 0.62F, 0.83F, 0.72F, 0.96F, 0.99F, 0.4F,  0.79F, 0.83F,
    0.21F, 0.43F, 0.32F, 0.44F, 0.72F, 0.21F, 0.4F,  0.93F, 0.07F, 0.11F, 0.41F, 0.24F, 0.04F,
    0.36F, 0.15F, 0.92F, 0.08F, 0.99F, 0.35F, 0.42F, 0.7F,  0.17F, 0.39F, 0.69F, 0.0F,  0.86F,
    0.89F, 0.59F, 0.81F, 0.77F, 0.15F, 0.89F, 0.17F, 0.76F, 0.67F, 0.58F, 0.78F, 0.26F, 0.19F,
    0.69F, 0.18F, 0.46F, 0.6F,  0.69F, 0.23F, 0.34F, 0.3F,  0.64F, 0.34F, 0.6F,  0.99F, 0.69F,
    0.57F, 0.75F, 0.07F, 0.36F, 0.75F, 0.81F, 0.8F,  0.42F, 0.09F, 0.94F, 0.66F, 0.35F, 0.67F,
    0.34F, 0.66F, 0.02F, 0.47F, 0.78F, 0.21F, 0.02F, 0.18F, 0.42F, 0.2F,  0.46F, 0.34F, 0.4F,
    0.46F, 0.96F, 0.86F, 0.25F, 0.25F, 0.22F, 0.37F, 0.59F, 0.19F, 0.45F, 0.61F, 0.04F, 0.71F,
    0.77F, 0.51F, 0.77F, 0.15F, 0.78F, 0.36F, 0.62F, 0.24F, 0.86F, 0.2F,  0.77F, 0.08F, 0.09F,
    0.3F,  0.0F,  0.6F,  0.99F, 0.69F,
};

//! A lot of bogus input spline derivatives - should have at list (max PME order = 5) * (DIM = 3) * (total unique atom number in all test cases = 16) values
std::vector<real> const c_sampleSplineDerivativesFull{
    0.82F, 0.88F, 0.83F, 0.11F, 0.93F, 0.32F, 0.71F, 0.37F, 0.69F, 0.88F, 0.11F, 0.38F, 0.25F,
    0.5F,  0.36F, 0.81F, 0.78F, 0.31F, 0.66F, 0.32F, 0.27F, 0.35F, 0.53F, 0.83F, 0.08F, 0.08F,
    0.94F, 0.71F, 0.65F, 0.24F, 0.13F, 0.01F, 0.33F, 0.65F, 0.24F, 0.53F, 0.45F, 0.84F, 0.33F,
    0.97F, 0.31F, 0.7F,  0.03F, 0.31F, 0.41F, 0.76F, 0.12F, 0.3F,  0.57F, 0.65F, 0.87F, 0.99F,
    0.42F, 0.97F, 0.32F, 0.39F, 0.73F, 0.23F, 0.03F, 0.67F, 0.97F, 0.57F, 0.42F, 0.38F, 0.54F,
    0.17F, 0.53F, 0.54F, 0.18F, 0.8F,  0.76F, 0.13F, 0.29F, 0.83F, 0.77F, 0.56F, 0.4F,  0.87F,
    0.36F, 0.18F, 0.59F, 0.04F, 0.05F, 0.61F, 0.26F, 0.91F, 0.62F, 0.16F, 0.89F, 0.23F, 0.26F,
    0.59F, 0.33F, 0.2F,  0.49F, 0.41F, 0.25F, 0.4F,  0.16F, 0.83F, 0.44F, 0.82F, 0.21F, 0.95F,
    0.14F, 0.8F,  0.37F, 0.31F, 0.41F, 0.53F, 0.15F, 0.85F, 0.78F, 0.17F, 0.92F, 0.03F, 0.13F,
    0.2F,  0.03F, 0.33F, 0.87F, 0.38F, 0,     0.08F, 0.79F, 0.36F, 0.53F, 0.05F, 0.07F, 0.94F,
    0.23F, 0.85F, 0.13F, 0.27F, 0.23F, 0.22F, 0.26F, 0.38F, 0.15F, 0.48F, 0.18F, 0.33F, 0.23F,
    0.62F, 0.1F,  0.36F, 0.99F, 0.07F, 0.02F, 0.04F, 0.09F, 0.29F, 0.52F, 0.29F, 0.83F, 0.97F,
    0.61F, 0.81F, 0.49F, 0.56F, 0.08F, 0.09F, 0.03F, 0.65F, 0.46F, 0.1F,  0.06F, 0.06F, 0.39F,
    0.29F, 0.04F, 0.03F, 0.1F,  0.83F, 0.94F, 0.59F, 0.97F, 0.82F, 0.2F,  0.66F, 0.23F, 0.11F,
    0.03F, 0.16F, 0.27F, 0.53F, 0.94F, 0.46F, 0.43F, 0.29F, 0.97F, 0.64F, 0.46F, 0.37F, 0.43F,
    0.48F, 0.37F, 0.93F, 0.5F,  0.2F,  0.92F, 0.09F, 0.74F, 0.55F, 0.44F, 0.05F, 0.13F, 0.17F,
    0.79F, 0.44F, 0.11F, 0.6F,  0.64F, 0.05F, 0.96F, 0.3F,  0.45F, 0.47F, 0.42F, 0.74F, 0.91F,
    0.06F, 0.89F, 0.24F, 0.26F, 0.68F, 0.4F,  0.88F, 0.5F,  0.65F, 0.48F, 0.15F, 0.0F,  0.41F,
    0.67F, 0.4F,  0.31F, 0.73F, 0.77F, 0.36F, 0.26F, 0.74F, 0.46F, 0.56F, 0.78F, 0.92F, 0.32F,
    0.9F,  0.06F, 0.55F, 0.6F,  0.13F, 0.38F, 0.93F, 0.5F,  0.92F, 0.96F, 0.82F, 0.0F,  0.04F,
    0.9F,  0.55F, 0.97F, 1.0F,  0.23F, 0.46F, 0.52F, 0.49F, 0.0F,  0.32F, 0.16F, 0.4F,  0.62F,
    0.36F, 0.03F, 0.63F, 0.16F, 0.58F, 0.97F, 0.03F, 0.44F, 0.07F, 0.22F, 0.75F, 0.32F, 0.61F,
    0.94F, 0.33F, 0.7F,  0.57F, 0.5F,  0.84F, 0.7F,  0.47F, 0.18F, 0.09F, 0.25F, 0.77F, 0.94F,
    0.85F, 0.09F, 0.83F, 0.02F, 0.91F,
};

//! 2 c_sample grids - only non-zero values have to be listed
std::vector<SparseRealGridValuesInput> const c_sampleGrids{ SparseRealGridValuesInput{
                                                                    { IVec{ 0, 0, 0 }, 3.5F },
                                                                    { IVec{ 7, 0, 0 }, -2.5F },
                                                                    { IVec{ 3, 5, 7 }, -0.006F },
                                                                    { IVec{ 1, 6, 7 }, -2.76F },
                                                                    { IVec{ 3, 1, 2 }, 0.6F },
                                                                    { IVec{ 6, 2, 4 }, 7.1F },
                                                                    { IVec{ 4, 5, 6 }, 4.1F },
                                                                    { IVec{ 4, 4, 6 }, -3.7F },
                                                            },
                                                            SparseRealGridValuesInput{
                                                                    { IVec{ 0, 4, 0 }, 6.F },
                                                                    { IVec{ 4, 2, 7 }, 13.76F },
                                                                    { IVec{ 0, 6, 7 }, 3.6F },
                                                                    { IVec{ 3, 2, 8 }, 0.61F },
                                                                    { IVec{ 5, 4, 3 }, 2.1F },
                                                                    { IVec{ 2, 5, 10 }, 3.6F },
                                                                    { IVec{ 5, 3, 6 }, 2.1F },
                                                                    { IVec{ 6, 6, 6 }, 6.1F },
                                                            } };

//! Input forces for reduction
std::vector<RVec> const c_sampleForcesFull{
    RVec{ 0.02F, 0.87F, 0.95F }, RVec{ 0.66F, 0.67F, 0.38F }, RVec{ 0.45F, 0.04F, 0.94F },
    RVec{ 0.54F, 0.76F, 0.58F }, RVec{ 0.83F, 0.31F, 0.73F }, RVec{ 0.71F, 0.06F, 0.35F },
    RVec{ 0.32F, 0.35F, 0.61F }, RVec{ 0.27F, 0.98F, 0.83F }, RVec{ 0.11F, 0.3F, 0.42F },
    RVec{ 0.95F, 0.69F, 0.58F }, RVec{ 0.29F, 0.1F, 0.68F },  RVec{ 0.94F, 0.62F, 0.51F },
    RVec{ 0.47F, 0.04F, 0.47F }, RVec{ 0.34F, 0.71F, 0.52F }
};

//! PME orders to test
std::vector<int> const pmeOrders{ 3, 4, 5 };
//! Atom counts to test
std::vector<size_t> const atomCounts{ 1, 2, 13 };

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
    GridLineIndicesVector gridLineIndices;
    //! Charges
    ChargesVector charges;
    //! Coordinates
    CoordinatesVector coordinates;
    //! Spline data for different orders
    std::map<int, AtomAndPmeOrderSizedData> splineDataByPmeOrder;
};

//! A range of the test input data sets, uniquely identified by the atom count
typedef std::map<size_t, AtomSizedData> InputDataByAtomCount;

/*! \brief Convenience typedef of the test input parameters - unit cell box, PME interpolation
 * order, grid dimensions, grid values, overwriting/reducing the input forces, atom count.
 *
 * The rest of the atom-related input data - gridline indices, spline theta values, spline dtheta
 * values, atom charges - is looked up in the inputAtomDataSets_ test fixture variable.
 */
typedef std::tuple<Matrix3x3, int, IVec, SparseRealGridValuesInput, PmeForceOutputHandling, size_t> GatherInputParameters;

//! Test fixture
class PmeGatherTest : public ::testing::TestWithParam<GatherInputParameters>
{
private:
    //! Storage of all the input atom datasets
    static InputDataByAtomCount s_inputAtomDataSets_;

public:
    PmeGatherTest() = default;
    //! Sets the input atom data references once
    static void SetUpTestCase()
    {
        size_t start = 0;
        for (auto atomCount : atomCounts)
        {
            AtomSizedData atomData;
            atomData.gridLineIndices =
                    GridLineIndicesVector(c_sampleGridLineIndicesFull).subArray(start, atomCount);
            atomData.charges = ChargesVector(c_sampleChargesFull).subArray(start, atomCount);
            start += atomCount;
            atomData.coordinates.resize(atomCount, RVec{ 1e6, 1e7, -1e8 });
            /* The coordinates are intentionally bogus in this test - only the size matters; the gridline indices are fed directly as inputs */
            for (auto pmeOrder : pmeOrders)
            {
                AtomAndPmeOrderSizedData splineData;
                const size_t             dimSize = atomCount * pmeOrder;
                for (int dimIndex = 0; dimIndex < DIM; dimIndex++)
                {
                    splineData.splineValues[dimIndex] =
                            SplineParamsDimVector(c_sampleSplineValuesFull).subArray(dimIndex * dimSize, dimSize);
                    splineData.splineDerivatives[dimIndex] =
                            SplineParamsDimVector(c_sampleSplineDerivativesFull).subArray(dimIndex * dimSize, dimSize);
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

        TestReferenceData refData;
        for (const auto& context : getPmeTestEnv()->getHardwareContexts())
        {
            CodePath   codePath = context->getCodePath();
            const bool supportedInput =
                    pmeSupportsInputForMode(*getPmeTestEnv()->hwinfo(), &inputRec, codePath);
            if (!supportedInput)
            {
                /* Testing the failure for the unsupported input */
                EXPECT_THROW(pmeInitWrapper(&inputRec, codePath, nullptr, nullptr, box), NotImplementedError);
                continue;
            }

            /* Describing the test uniquely */
            SCOPED_TRACE(
                    formatString("Testing force gathering with %s %sfor PME grid size %d %d %d"
                                 ", order %d, %zu atoms, %s",
                                 codePathToString(codePath), context->getDescription().c_str(),
                                 gridSize[XX], gridSize[YY], gridSize[ZZ], pmeOrder, atomCount,
                                 (inputForceTreatment == PmeForceOutputHandling::ReduceWithInput)
                                         ? "with reduction"
                                         : "without reduction"));

            PmeSafePointer pmeSafe = pmeInitWrapper(&inputRec, codePath, context->getDeviceInfo(),
                                                    context->getPmeGpuProgram(), box);
            std::unique_ptr<StatePropagatorDataGpu> stateGpu =
                    (codePath == CodePath::GPU) ? makeStatePropagatorDataGpu(*pmeSafe.get()) : nullptr;

            pmeInitAtoms(pmeSafe.get(), stateGpu.get(), codePath, inputAtomData.coordinates,
                         inputAtomData.charges);

            /* Setting some more inputs */
            pmeSetRealGrid(pmeSafe.get(), codePath, nonZeroGridValues);
            pmeSetGridLineIndices(pmeSafe.get(), codePath, inputAtomData.gridLineIndices);
            for (int dimIndex = 0; dimIndex < DIM; dimIndex++)
            {
                pmeSetSplineData(pmeSafe.get(), codePath, inputAtomSplineData.splineValues[dimIndex],
                                 PmeSplineDataType::Values, dimIndex);
                pmeSetSplineData(pmeSafe.get(), codePath, inputAtomSplineData.splineDerivatives[dimIndex],
                                 PmeSplineDataType::Derivatives, dimIndex);
            }

            /* Explicitly copying the sample forces to be able to modify them */
            auto inputForcesFull(c_sampleForcesFull);
            GMX_RELEASE_ASSERT(inputForcesFull.size() >= atomCount, "Bad input forces size");
            auto forces = ForcesVector(inputForcesFull).subArray(0, atomCount);

            /* Running the force gathering itself */
            pmePerformGather(pmeSafe.get(), codePath, inputForceTreatment, forces);
            pmeFinalizeTest(pmeSafe.get(), codePath);

            /* Check the output forces correctness */
            TestReferenceChecker forceChecker(refData.rootChecker());
            const auto           ulpTolerance = 3 * pmeOrder;
            forceChecker.setDefaultTolerance(relativeToleranceAsUlp(1.0, ulpTolerance));
            forceChecker.checkSequence(forces.begin(), forces.end(), "Forces");
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
INSTANTIATE_TEST_CASE_P(SaneInput,
                        PmeGatherTest,
                        ::testing::Combine(::testing::ValuesIn(c_sampleBoxes),
                                           ::testing::ValuesIn(pmeOrders),
                                           ::testing::ValuesIn(c_sampleGridSizes),
                                           ::testing::ValuesIn(c_sampleGrids),
                                           ::testing::Values(PmeForceOutputHandling::Set,
                                                             PmeForceOutputHandling::ReduceWithInput),
                                           ::testing::ValuesIn(atomCounts)));

} // namespace
} // namespace test
} // namespace gmx
