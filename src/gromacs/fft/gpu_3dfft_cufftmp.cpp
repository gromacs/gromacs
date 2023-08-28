/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 *  \brief Implements GPU 3D FFT routines using HeFFTe.
 *
 *  \author Gaurav Garg <gaugarg@nvidia.com>
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_cufftmp.h"

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#if CUFFT_VERSION >= 11005
// Extra interface which maps legacy cufftXtSetDistribution API to the updated version
// first introduced in the NVIDIA HPC SDK 23.3 (cuFFT/cuFFTMp version 11.0.5).

/*! \brief Definition of box type for input or output of FFT operation, with lower and upper indices and strides */
typedef struct cufftBox3d_t
{
    //! Lower indices in each of 3 dimensions for 3D box
    size_t lower[3];
    //! Upper indices in each of 3 dimensions for 3D box
    size_t upper[3];
    //! Strides in each of 3 dimensions for 3D box
    size_t strides[3];
} cufftBox3d;

/*! \brief
 * Interface the legacy cufftXtSetDistribution call to the new version
 * \param [in] plan         Plan for FFT
 * \param [in] inputBox     Input box for FFT
 * \param [in] outputBox    Output box for FFT
 * \returns cufftResult_t   Result code
 */
static cufftResult_t cufftXtSetDistribution(cufftHandle plan, const cufftBox3d* inputBox, cufftBox3d* outputBox)
{

    long long int inputLower[3];
    long long int inputUpper[3];
    long long int inputStrides[3];
    long long int outputLower[3];
    long long int outputUpper[3];
    long long int outputStrides[3];

    for (int i = 0; i < 3; i++)
    {
        inputLower[i]    = inputBox->lower[i];
        inputUpper[i]    = inputBox->upper[i];
        inputStrides[i]  = inputBox->strides[i];
        outputLower[i]   = outputBox->lower[i];
        outputUpper[i]   = outputBox->upper[i];
        outputStrides[i] = outputBox->strides[i];
    }

    cufftResult_t result = cufftXtSetDistribution(
            plan, 3, inputLower, inputUpper, outputLower, outputUpper, inputStrides, outputStrides);

    return result;
}
#endif

namespace gmx
{
static void handleCufftError(cufftResult_t status, const char* msg)
{
    if (status != CUFFT_SUCCESS)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("%s (error code %d)\n", msg, status)));
    }
}

Gpu3dFft::ImplCuFftMp::ImplCuFftMp(bool                allocateRealGrid,
                                   MPI_Comm            comm,
                                   ArrayRef<const int> gridSizesInXForEachRank,
                                   ArrayRef<const int> gridSizesInYForEachRank,
                                   const int           nz,
                                   bool                performOutOfPlaceFFT,
                                   const DeviceContext& /*context*/,
                                   const DeviceStream&  pmeStream,
                                   ivec                 realGridSize,
                                   ivec                 realGridSizePadded,
                                   ivec                 complexGridSizePadded,
                                   DeviceBuffer<float>* realGrid,
                                   DeviceBuffer<float>* complexGrid)
{
    const int numDomainsX = gridSizesInXForEachRank.size();
    const int numDomainsY = gridSizesInYForEachRank.size();

    GMX_RELEASE_ASSERT(allocateRealGrid == true, "Grids cannot be pre-allocated");
    GMX_RELEASE_ASSERT(performOutOfPlaceFFT == false, "Only in-place FFT supported");
    GMX_RELEASE_ASSERT(numDomainsX * numDomainsY > 1,
                       "CuFftMp backend is expected to be used only with more than 1 rank");

    // calculate grid offsets
    std::vector<size_t> gridOffsetsInX(numDomainsX + 1);
    std::vector<size_t> gridOffsetsInY(numDomainsY + 1);

    gridOffsetsInX[0] = 0;
    for (unsigned int i = 0; i < gridSizesInXForEachRank.size(); ++i)
    {
        gridOffsetsInX[i + 1] = gridOffsetsInX[i] + gridSizesInXForEachRank[i];
    }

    gridOffsetsInY[0] = 0;
    for (unsigned int i = 0; i < gridSizesInYForEachRank.size(); ++i)
    {
        gridOffsetsInY[i + 1] = gridOffsetsInY[i] + gridSizesInYForEachRank[i];
    }

    int rank, nProcs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nProcs);

    GMX_RELEASE_ASSERT(nProcs == numDomainsX * numDomainsY,
                       "Mismatch in communicator size and expected domain decomposition");

    // define how ranks are mapped to 2d domain
    int procY = rank % numDomainsY;
    int procX = rank / numDomainsY;

    // output dimension padded along z
    const size_t complexZDim = nz / 2 + 1;
    const size_t nzPadded    = complexZDim * 2;

    // local real grid boxes - pencil in X & Y, along Z
    cufftBox3d const realBox = { { gridOffsetsInX[procX], gridOffsetsInY[procY], 0 },
                                 { gridOffsetsInX[procX + 1], gridOffsetsInY[procY + 1], (size_t)nz },
                                 { gridSizesInYForEachRank[procY] * nzPadded, nzPadded, 1 } };

    const size_t nx = gridOffsetsInX[numDomainsX];
    const size_t ny = gridOffsetsInY[numDomainsY];

    cufftBox3d complexBox;
    size_t     gridSizesInY_transformed = ny;
    size_t     gridSizesInZ_transformed = complexZDim;

    // if possible, keep complex data in slab decomposition along Z
    // this allows cuFFTMp to have single communication phase
    if (numDomainsY > 1 && complexZDim >= (size_t)nProcs)
    {
        // define shape of local complex grid boxes
        std::vector<size_t> gridOffsetsInZ_transformed(nProcs + 1);

        for (int i = 0; i < nProcs; i++)
        {
            gridOffsetsInZ_transformed[i] = (i * complexZDim) / nProcs;
        }
        gridOffsetsInZ_transformed[nProcs] = complexZDim;

        gridSizesInZ_transformed =
                gridOffsetsInZ_transformed[rank + 1] - gridOffsetsInZ_transformed[rank];

        // Output data are slab decomposed along Z
        complexBox = { { 0, 0, gridOffsetsInZ_transformed[rank] },
                       { nx, ny, gridOffsetsInZ_transformed[rank + 1] },
                       { ny * gridSizesInZ_transformed, gridSizesInZ_transformed, 1 } };
    }
    else
    {
        // define shape of local complex grid boxes
        std::vector<size_t> gridOffsetsInY_transformed(numDomainsX + 1);
        std::vector<size_t> gridOffsetsInZ_transformed(numDomainsY + 1);

        for (int i = 0; i < numDomainsX; i++)
        {
            gridOffsetsInY_transformed[i] = (i * ny + 0) / numDomainsX;
        }
        gridOffsetsInY_transformed[numDomainsX] = ny;

        for (int i = 0; i < numDomainsY; i++)
        {
            gridOffsetsInZ_transformed[i] = (i * complexZDim + 0) / numDomainsY;
        }
        gridOffsetsInZ_transformed[numDomainsY] = complexZDim;

        gridSizesInY_transformed =
                gridOffsetsInY_transformed[procX + 1] - gridOffsetsInY_transformed[procX];
        gridSizesInZ_transformed =
                gridOffsetsInZ_transformed[procY + 1] - gridOffsetsInZ_transformed[procY];

        // Output data are pencils in Y & Z, along X
        complexBox = {
            { 0, gridOffsetsInY_transformed[procX], gridOffsetsInZ_transformed[procY] },
            { nx, gridOffsetsInY_transformed[procX + 1], gridOffsetsInZ_transformed[procY + 1] },
            { gridSizesInY_transformed * gridSizesInZ_transformed, gridSizesInZ_transformed, 1 }
        };
    }

    handleCufftError(cufftCreate(&planR2C_), "cufftCreate R2C plan failure");
    handleCufftError(cufftCreate(&planC2R_), "cufftCreate C2R plan failure");

    // Duplicate the existing communicator for cuFFTMp usage as the communicator should
    // remain valid from cuFFTMp plan creation to destruction and should remain in scope.
    MPI_Comm_dup(comm, &comm_);

    // Attach the MPI communicator to the plans
    handleCufftError(cufftMpAttachComm(planR2C_, CUFFT_COMM_MPI, &comm_),
                     "cufftMpAttachComm R2C plan failure");
    handleCufftError(cufftMpAttachComm(planC2R_, CUFFT_COMM_MPI, &comm_),
                     "cufftMpAttachComm C2R plan failure");

    // Describe the data distribution
    handleCufftError(cufftXtSetDistribution(planR2C_, &realBox, &complexBox),
                     "cufftXtSetDistribution R2C plan failure");
    handleCufftError(cufftXtSetDistribution(planC2R_, &realBox, &complexBox),
                     "cufftXtSetDistribution C2R plan failure");

    // Set the stream
    handleCufftError(cufftSetStream(planR2C_, pmeStream.stream()),
                     "cufftSetStream R2C plan failure");
    handleCufftError(cufftSetStream(planC2R_, pmeStream.stream()),
                     "cufftSetStream C2R plan failure");

    // Make the plan
    size_t workspace;
    handleCufftError(cufftMakePlan3d(planR2C_, nx, ny, nz, CUFFT_R2C, &workspace),
                     "cufftMakePlan3d R2C plan failure");
    handleCufftError(cufftMakePlan3d(planC2R_, nx, ny, nz, CUFFT_C2R, &workspace),
                     "cufftMakePlan3d C2R plan failure");

    // Allocate GPU memory
    // Data is initially distributed according to CUFFT_XT_FORMAT_INPLACE
    handleCufftError(cufftXtMalloc(planR2C_, &desc_, CUFFT_XT_FORMAT_DISTRIBUTED_INPUT),
                     "cufftXtMalloc failure");

    // write back the output data
    *realGrid    = (float*)desc_->descriptor->data[0];
    *complexGrid = (float*)desc_->descriptor->data[0];

    realGridSize[XX] = gridSizesInXForEachRank[procX];
    realGridSize[YY] = gridSizesInYForEachRank[procY];
    realGridSize[ZZ] = nz;

    realGridSizePadded[XX] = gridSizesInXForEachRank[procX];
    realGridSizePadded[YY] = gridSizesInYForEachRank[procY];
    realGridSizePadded[ZZ] = nzPadded;

    complexGridSizePadded[XX] = nx;
    complexGridSizePadded[YY] = gridSizesInY_transformed;
    complexGridSizePadded[ZZ] = gridSizesInZ_transformed;
}

Gpu3dFft::ImplCuFftMp::~ImplCuFftMp()
{
    handleCufftError(cufftXtFree(desc_), "cufftXtFree failure");

    handleCufftError(cufftDestroy(planR2C_), "cufftDestroy R2C plan failure");
    handleCufftError(cufftDestroy(planC2R_), "cufftDestroy C2R plan failure");
    MPI_Comm_free(&comm_);
}

void Gpu3dFft::ImplCuFftMp::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
    switch (dir)
    {
        case GMX_FFT_REAL_TO_COMPLEX:
            // Run R2C
            handleCufftError(cufftXtExecDescriptor(planR2C_, desc_, desc_, CUFFT_FORWARD),
                             "cufftXtExecDescriptor R2C plan failure");
            break;
        case GMX_FFT_COMPLEX_TO_REAL:
            // Run C2R
            handleCufftError(cufftXtExecDescriptor(planC2R_, desc_, desc_, CUFFT_INVERSE),
                             "cufftXtExecDescriptor C2R plan failure");
            break;
        default:
            GMX_THROW(NotImplementedError("The chosen 3D-FFT case is not implemented on GPUs"));
    }
}
} // namespace gmx
