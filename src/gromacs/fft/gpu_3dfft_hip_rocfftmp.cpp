/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 *  \brief Implements GPU 3D FFT routines for HIP via rocFFT.
 *
 *  This variant includes MPI distribution of work to different ranks
 *  and GPU devices.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_hip_rocfftmp.h"

#include <array>
#include <type_traits>
#include <vector>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fixedcapacityvector.h"
#include "gromacs/utility/gmxassert.h"

#include "rocfft_common_utils.h"

#if !defined(__HIPCC__)
#    error This file can only be compiled with the HIPCC compiler for ROCm
#endif

namespace gmx
{

namespace
{

//! Number of dimensions plus one for batch size.
static constexpr int sc_dimPlusBatch = DIM + 1;

/*! \brief Definition of box type for input or output of FFT operation.
 *
 * Maps to data required by the rocfft_brick type to setup the decomposition.
 * with lower and upper indices and strides */
struct rocFFTBoxType
{
    //! Lower indices in each dimensions
    size_t lower[sc_dimPlusBatch];
    //! Upper indices in each dimension
    size_t upper[sc_dimPlusBatch];
    //! Strides in each dimensions
    size_t strides[sc_dimPlusBatch];
    //! Get size of box in \p dim
    size_t size(int dim) const
    {
        GMX_RELEASE_ASSERT(dim < DIM, "index out of range");
        return upper[dim] - lower[dim];
    }
};

size_t calculateBoxSize(const rocFFTBoxType& box)
{
    // Calculate buffer size based on the stride-aware memory layout.
    // The box dimensions are ordered as [Z, Y, X, batch] (rocFFT convention)
    // but strides encode GROMACS [X][Y][Z] memory layout (Z fastest-varying).
    // The total size is the extent of the slowest-varying dimension (X, index 2)
    // multiplied by its stride, which accounts for all padding in faster dimensions.
    constexpr int xDimension = 2; // X is at index 2 in [Z, Y, X, batch] ordering
    const size_t  extentInX  = box.upper[xDimension] - box.lower[xDimension];
    const size_t  strideInX  = box.strides[xDimension];
    return extentInX * strideInX;
}

/*! \brief Prepare plans for the forward and reverse transformation.
 *
 * Because these require device-side allocations, some of them must be
 * done in a SYCL queue. */
RocfftPlan makePlan(const char*            descriptiveString,
                    rocfft_transform_type  transformType,
                    rocfft_array_type      inputArrayType,
                    rocfft_array_type      outputArrayType,
                    ArrayRef<const size_t> rocfftRealGridSize,
                    MPI_Comm               comm,
                    int                    deviceId,
                    rocFFTBoxType          inputGridBox,
                    rocFFTBoxType          outputGridBox)
{
    rocfft_plan_description description = nullptr;
    rocfft_status           result;
    result = rocfft_plan_description_create(&description);
    handleRocFftError(result, descriptiveString, "rocfft_plan_description_create failure");
    result = rocfft_plan_description_set_comm(description, rocfft_comm_mpi, &comm);
    handleRocFftError(result, descriptiveString, "rocfft_plan_description_set_comm failure");
    rocfft_field inputField;
    result = rocfft_field_create(&inputField);
    handleRocFftError(result, descriptiveString, "rocfft_field_create failure for input field");
    rocfft_field outputField;
    result = rocfft_field_create(&outputField);
    handleRocFftError(result, descriptiveString, "rocfft_field_create failure for output field");


    rocfft_brick inputGridBrick;
    result = rocfft_brick_create(
            &inputGridBrick, inputGridBox.lower, inputGridBox.upper, inputGridBox.strides, sc_dimPlusBatch, deviceId);
    handleRocFftError(result, descriptiveString, "rocfft_brick_create failure for inputGridBrick");
    result = rocfft_field_add_brick(inputField, inputGridBrick);
    handleRocFftError(result,
                      descriptiveString,
                      "rocfft_field_add_brick failure for inputGridBrick to inputField");
    result = rocfft_brick_destroy(inputGridBrick);
    handleRocFftError(result, descriptiveString, "rocfft_brick_destroy failure for inputGridBrick");

    rocfft_brick outputGridBrick;
    result = rocfft_brick_create(&outputGridBrick,
                                 outputGridBox.lower,
                                 outputGridBox.upper,
                                 outputGridBox.strides,
                                 sc_dimPlusBatch,
                                 deviceId);
    handleRocFftError(result, descriptiveString, "rocfft_brick_create failure for outputGridBrick");
    result = rocfft_field_add_brick(outputField, outputGridBrick);
    handleRocFftError(result,
                      descriptiveString,
                      "rocfft_field_add_brick failure for outputGridBrick to outputField");
    result = rocfft_brick_destroy(outputGridBrick);
    handleRocFftError(result, descriptiveString, "rocfft_brick_destroy failure for outputGridBrick");

    result = rocfft_plan_description_add_infield(description, inputField);
    handleRocFftError(result, descriptiveString, "rocfft_plan_description_add_infield failure");
    result = rocfft_field_destroy(inputField);
    handleRocFftError(result, descriptiveString, "rocfft_field_destroy failure for inputField");

    result = rocfft_plan_description_add_outfield(description, outputField);
    handleRocFftError(result, descriptiveString, "rocfft_plan_description_add_outfield failure");
    result = rocfft_field_destroy(outputField);
    handleRocFftError(result, descriptiveString, "rocfft_field_destroy failure for outputField");

    result = rocfft_plan_description_set_data_layout(description,
                                                     inputArrayType,
                                                     outputArrayType,
                                                     // No offsets are needed
                                                     nullptr,
                                                     nullptr,
                                                     0,
                                                     nullptr,
                                                     0,
                                                     0,
                                                     nullptr,
                                                     0);
    handleRocFftError(result, descriptiveString, "rocfft_plan_description_set_data_layout failure");

    size_t requiredWorkBufferSize = 0;
    void*  workBuffer             = nullptr;

    // The plan creation depends on the identity of the GPU device, so
    // we make sure it is made in the same queue where it will be
    // used. The stream for execution can be set at the same time.

    // First set up device buffers to receive the rocfft status values
    rocfft_plan   plan;
    constexpr int numBatches = 1;
    result                   = rocfft_plan_create(&plan,
                                rocfft_placement_notinplace,
                                transformType,
                                rocfft_precision_single,
                                rocfftRealGridSize.size(),
                                rocfftRealGridSize.data(),
                                numBatches,
                                description);
    handleRocFftError(result, descriptiveString, "rocfft_plan_create failure");
    result = rocfft_plan_description_destroy(description);
    handleRocFftError(result, descriptiveString, "rocfft_plan_description_destroy failure");

    result = rocfft_plan_get_work_buffer_size(plan, &requiredWorkBufferSize);
    handleRocFftError(result, descriptiveString, "rocfft_plan_get_work_buffer_size failure");
    if (requiredWorkBufferSize > 0)
    {
        hipError_t err = hipMalloc(&workBuffer, requiredWorkBufferSize);
        result         = (err == hipSuccess) ? rocfft_status_success : rocfft_status_failure;
    }
    else
    {
        result = rocfft_status_success;
    }
    handleRocFftError(result, descriptiveString, "hipMalloc failure");

    rocfft_execution_info execution_info = nullptr;
    result                               = rocfft_execution_info_create(&execution_info);
    handleRocFftError(result, descriptiveString, "rocfft_execution_info_create failure");

    if (requiredWorkBufferSize > 0)
    {
        GMX_RELEASE_ASSERT(workBuffer != nullptr,
                           "Work buffer should have been allocated, but was not");
        result = rocfft_execution_info_set_work_buffer(execution_info, workBuffer, requiredWorkBufferSize);
        handleRocFftError(
                result, descriptiveString, "rocfft_execution_info_set_work_buffer failure");
    }

    return RocfftPlan{ plan, execution_info, workBuffer };
}

} // namespace

//! Impl class
class Gpu3dFft::ImplHipRocfftMp::Impl
{
public:
    //! \copydoc Gpu3dFft::Impl::Impl
    Impl(MPI_Comm            comm,
         int                 deviceId,
         const DeviceStream& pmeStream,
         ArrayRef<size_t>    realGridSize,
         rocFFTBoxType       realBox,
         rocFFTBoxType       complexBox);

    ~Impl();

    /*! \brief Handle initializing the rocFFT library
     *
     * Make sure the library is initialized before the plans, etc. and
     * not destructed before they are. */
    RocfftInitializer init_;
    //! Data for 3D FFT plans and execution
    EnumerationArray<FftDirection, RocfftPlan> plans_;
    const DeviceStream&                        pmeStream_;
};

Gpu3dFft::ImplHipRocfftMp::Impl::Impl(MPI_Comm            comm,
                                      const int           deviceId,
                                      const DeviceStream& pmeStream,
                                      ArrayRef<size_t>    realGridSize,
                                      rocFFTBoxType       realBox,
                                      rocFFTBoxType       complexBox) :
    plans_({ makePlan("real-to-complex",
                      rocfft_transform_type_real_forward,
                      rocfft_array_type_real,
                      rocfft_array_type_hermitian_interleaved,
                      realGridSize,
                      comm,
                      deviceId,
                      realBox,
                      complexBox),
             makePlan("complex-to-real",
                      rocfft_transform_type_real_inverse,
                      rocfft_array_type_hermitian_interleaved,
                      rocfft_array_type_real,
                      realGridSize,
                      comm,
                      deviceId,
                      complexBox,
                      realBox) }),
    pmeStream_(pmeStream)
{
}

Gpu3dFft::ImplHipRocfftMp::Impl::~Impl() {}

void Gpu3dFft::ImplHipRocfftMp::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
    GMX_RELEASE_ASSERT(
            (dir == GMX_FFT_REAL_TO_COMPLEX) || (dir == GMX_FFT_COMPLEX_TO_REAL),
            "Only real-to-complex and complex-to-real FFTs are implemented in HIP rocfft");
    FftDirection         direction;
    std::array<void*, 1> inputGridPointers, outputGridPointers;

    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        direction             = FftDirection::RealToComplex;
        inputGridPointers[0]  = reinterpret_cast<void*>(realGrid_);
        outputGridPointers[0] = reinterpret_cast<void*>(complexGrid_);
    }
    else
    {
        direction             = FftDirection::ComplexToReal;
        inputGridPointers[0]  = reinterpret_cast<void*>(complexGrid_);
        outputGridPointers[0] = reinterpret_cast<void*>(realGrid_);
    }
    rocfft_execution_info_set_stream(impl_->plans_[direction].info, impl_->pmeStream_.stream());
    // Don't check results generated asynchronously,
    // because we don't know what to do with them
    rocfft_execute(impl_->plans_[direction].plan,
                   inputGridPointers.data(),
                   outputGridPointers.data(),
                   impl_->plans_[direction].info);
}

Gpu3dFft::ImplHipRocfftMp::ImplHipRocfftMp(bool                 allocateRealGrid,
                                           MPI_Comm             comm,
                                           ArrayRef<const int>  gridSizesInXForEachRank,
                                           ArrayRef<const int>  gridSizesInYForEachRank,
                                           const int            nz,
                                           bool                 performOutOfPlaceFFT,
                                           const DeviceContext& context,
                                           const DeviceStream&  pmeStream,
                                           ivec                 realGridSize,
                                           ivec                 realGridSizePadded,
                                           ivec                 complexGridSizePadded,
                                           DeviceBuffer<float>* realGrid,
                                           DeviceBuffer<float>* complexGrid) :
    Gpu3dFft::Impl::Impl(performOutOfPlaceFFT)
{
    MPI_Comm_dup(comm, &comm_);
    GMX_RELEASE_ASSERT(performOutOfPlaceFFT,
                       "Only out-of-place FFT is implemented in rocfft with decomposition");
    GMX_RELEASE_ASSERT(allocateRealGrid == true, "Grids need to be allocated here");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() * gridSizesInYForEachRank.size() > 1,
                       "Distributed rocfft expected to be used with more than 1 rank");

    const int numDomainsX = gridSizesInXForEachRank.size();
    const int numDomainsY = gridSizesInYForEachRank.size();

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

    // Use a packed real-grid layout. The in-place R2C Z-padding inherited
    // from the cuFFTMp adapter is wrong here: this plan is out-of-place
    // (see the GMX_RELEASE_ASSERT(performOutOfPlaceFFT, ...) above), and
    // the padded brick stride produced silently wrong PME values.
    const size_t complexZDim = nz / 2 + 1;
    const size_t nzReal      = static_cast<size_t>(nz);

    // local real grid boxes - pencil in X & Y, along Z
    const rocFFTBoxType realBox = { { 0, gridOffsetsInY[procY], gridOffsetsInX[procX], 0 },
                                    { nzReal, gridOffsetsInY[procY + 1], gridOffsetsInX[procX + 1], 1 },
                                    { 1, nzReal, gridSizesInYForEachRank[procY] * nzReal, 0 } };
    const size_t        nx      = gridOffsetsInX[numDomainsX];
    const size_t        ny      = gridOffsetsInY[numDomainsY];

    rocFFTBoxType complexBox;
    size_t        gridSizesInY_transformed = ny;
    size_t        gridSizesInZ_transformed = complexZDim;

    // if possible, keep complex data in slab decomposition along Z
    // this allows rocFFTMp to have single communication phase
    if (numDomainsY > 1 && complexZDim >= static_cast<size_t>(nProcs))
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

        complexBox = { { gridOffsetsInZ_transformed[rank], 0, 0, 0 },
                       { gridOffsetsInZ_transformed[rank + 1], ny, nx, 1 },
                       { 1, gridSizesInZ_transformed, ny * gridSizesInZ_transformed, 0 } };
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

        complexBox = {
            { gridOffsetsInZ_transformed[procY], gridOffsetsInY_transformed[procX], 0, 0 },
            { gridOffsetsInZ_transformed[procY + 1], gridOffsetsInY_transformed[procX + 1], nx, 1 },
            { 1, gridSizesInZ_transformed, gridSizesInY_transformed * gridSizesInZ_transformed, 0 }
        };
    }

    realGridSize[XX] = gridSizesInXForEachRank[procX];
    realGridSize[YY] = gridSizesInYForEachRank[procY];
    realGridSize[ZZ] = nz;

    realGridSizePadded[XX] = gridSizesInXForEachRank[procX];
    realGridSizePadded[YY] = gridSizesInYForEachRank[procY];
    realGridSizePadded[ZZ] = nz;

    complexGridSizePadded[XX] = nx;
    complexGridSizePadded[YY] = gridSizesInY_transformed;
    complexGridSizePadded[ZZ] = gridSizesInZ_transformed;

    std::vector<size_t> rocfftRealGridSize = { static_cast<size_t>(nz),
                                               static_cast<size_t>(ny),
                                               static_cast<size_t>(nx) };

    const size_t realSize    = calculateBoxSize(realBox);
    const size_t complexSize = calculateBoxSize(complexBox) * 2;

    allocateDeviceBuffer(&realGrid_, realSize, context);
    allocateDeviceBuffer(&complexGrid_, complexSize, context);

    *realGrid    = realGrid_;
    *complexGrid = complexGrid_;

    const int deviceId = context.deviceInfo().id;

    impl_ = std::make_unique<Impl>(comm_, deviceId, pmeStream, rocfftRealGridSize, realBox, complexBox);
}

Gpu3dFft::ImplHipRocfftMp::~ImplHipRocfftMp()
{
    MPI_Comm_free(&comm_);
    freeDeviceBuffer(&realGrid_);
    freeDeviceBuffer(&complexGrid_);
}

} // namespace gmx
