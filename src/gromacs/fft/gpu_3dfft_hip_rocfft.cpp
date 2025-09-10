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
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_hip_rocfft.h"

#include <vector>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "rocfft_common_utils.h"

#if !defined(__HIPCC__)
#    error This file can only be compiled with the HIPCC compiler for ROCm
#endif

namespace gmx
{

namespace
{

/*! \brief Prepare plans for the forward and reverse transformation.
 *
 * Because these require device-side allocations, some of them must be
 * done in a SYCL queue. */
RocfftPlan makePlan(const std::string&     descriptiveString,
                    rocfft_transform_type  transformType,
                    const PlanSetupData&   inputPlanSetupData,
                    const PlanSetupData&   outputPlanSetupData,
                    ArrayRef<const size_t> rocfftRealGridSize)
{
    rocfft_plan_description description = nullptr;
    rocfft_status           result;
    result = rocfft_plan_description_create(&description);
    handleRocFftError(result, descriptiveString, "rocfft_plan_description_create failure");
    result = rocfft_plan_description_set_data_layout(description,
                                                     inputPlanSetupData.arrayType,
                                                     outputPlanSetupData.arrayType,
                                                     // No offsets are needed
                                                     nullptr,
                                                     nullptr,
                                                     inputPlanSetupData.strides.size(),
                                                     inputPlanSetupData.strides.data(),
                                                     inputPlanSetupData.totalSize,
                                                     outputPlanSetupData.strides.size(),
                                                     outputPlanSetupData.strides.data(),
                                                     outputPlanSetupData.totalSize);
    handleRocFftError(result, descriptiveString, "rocfft_plan_description_set_data_layout failure");

    size_t requiredWorkBufferSize = 0;
    void*  workBuffer             = nullptr;


    // The plan creation depends on the identity of the GPU device, so
    // we make sure it is made in the same queue where it will be
    // used. The stream for execution can be set at the same time.

    // First set up device buffers to receive the rocfft status values
    rocfft_plan plan;
    const int   numBatches = 1;
    result                 = rocfft_plan_create(&plan,
                                rocfft_placement_notinplace,
                                transformType,
                                rocfft_precision_single,
                                rocfftRealGridSize.size(),
                                rocfftRealGridSize.data(),
                                numBatches,
                                description);
    handleRocFftError(result, descriptiveString, "rocfft_plan_create failure");

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

    result = rocfft_plan_description_destroy(description);
    handleRocFftError(result, descriptiveString, "rocfft_plan_description_destroy failure");

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
class Gpu3dFft::ImplHipRocfft::Impl
{
public:
    //! \copydoc Gpu3dFft::Impl::Impl
    Impl(bool                 allocateRealGrid,
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
         DeviceBuffer<float>* complexGrid);

    /*! \brief Handle initializing the rocFFT library
     *
     * Make sure the library is initialized before the plans, etc. and
     * not destructed before they are. */
    RocfftInitializer init_;
    //! Data for 3D FFT plans and execution
    EnumerationArray<FftDirection, RocfftPlan> plans_;
    //! Handle to the real grid buffer
    float*              realGrid_;
    float*              complexGrid_;
    const DeviceStream& pmeStream_;
};

Gpu3dFft::ImplHipRocfft::Impl::Impl(bool allocateRealGrid,
                                    MPI_Comm /*comm*/,
                                    ArrayRef<const int> gridSizesInXForEachRank,
                                    ArrayRef<const int> gridSizesInYForEachRank,
                                    int /*nz*/,
                                    bool performOutOfPlaceFFT,
                                    const DeviceContext& /*context*/,
                                    const DeviceStream&  pmeStream,
                                    ivec                 realGridSize,
                                    ivec                 realGridSizePadded,
                                    ivec                 complexGridSizePadded,
                                    DeviceBuffer<float>* realGrid,
                                    DeviceBuffer<float>* /*complexGrid*/) :
    plans_{ makePlan("real-to-complex",
                     rocfft_transform_type_real_forward,
                     // input
                     PlanSetupData{ rocfft_array_type_real,
                                    makeRealStrides(realGridSizePadded),
                                    computeTotalSize(realGridSizePadded) },
                     // output
                     PlanSetupData{ rocfft_array_type_hermitian_interleaved,
                                    makeComplexStrides(complexGridSizePadded),
                                    computeTotalSize(complexGridSizePadded) },
                     // Note that rocFFT requires that we reverse the dimension order when planning
                     std::vector<size_t>{ size_t(realGridSize[ZZ]),
                                          size_t(realGridSize[YY]),
                                          size_t(realGridSize[XX]) }),
            // For rocFFT, the complex-to-real setup is the logical
            // converse of the real-to-complex. The PlanSetupData objects
            // are the same, but used in the opposite sense of
            // input/output.
            makePlan("complex-to-real",
                     rocfft_transform_type_real_inverse,
                     // input
                     PlanSetupData{ rocfft_array_type_hermitian_interleaved,
                                    makeComplexStrides(complexGridSizePadded),
                                    computeTotalSize(complexGridSizePadded) },
                     // output
                     PlanSetupData{ rocfft_array_type_real,
                                    makeRealStrides(realGridSizePadded),
                                    computeTotalSize(realGridSizePadded) },
                     // Note that rocFFT requires that we reverse the dimension order when planning
                     std::vector<size_t>{ size_t(realGridSize[ZZ]),
                                          size_t(realGridSize[YY]),
                                          size_t(realGridSize[XX]) }) },
    realGrid_(*realGrid),
    pmeStream_(pmeStream)
{
    GMX_RELEASE_ASSERT(performOutOfPlaceFFT, "Only out-of-place FFT is implemented in rocfft");
    GMX_RELEASE_ASSERT(allocateRealGrid == false, "Grids need to be pre-allocated");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() == 1 && gridSizesInYForEachRank.size() == 1,
                       "FFT decomposition not implemented with the default HIP rocfft backend");
}

void Gpu3dFft::ImplHipRocfft::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
    GMX_RELEASE_ASSERT(
            (dir == GMX_FFT_REAL_TO_COMPLEX) || (dir == GMX_FFT_COMPLEX_TO_REAL),
            "Only real-to-complex and complex-to-real FFTs are implemented in HIP rocfft");
    FftDirection direction;
    float **     inputGrid = nullptr, **outputGrid = nullptr;
    impl_->complexGrid_ = complexGrid_;
    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        direction  = FftDirection::RealToComplex;
        inputGrid  = &impl_->realGrid_;
        outputGrid = &impl_->complexGrid_;
    }
    else
    {
        direction  = FftDirection::ComplexToReal;
        inputGrid  = &impl_->complexGrid_;
        outputGrid = &impl_->realGrid_;
    }
    void* d_inputGrid  = reinterpret_cast<void*>(*inputGrid);
    void* d_outputGrid = reinterpret_cast<void*>(*outputGrid);
    rocfft_execution_info_set_stream(impl_->plans_[direction].info, impl_->pmeStream_.stream());
    // Don't check results generated asynchronously,
    // because we don't know what to do with them
    rocfft_execute(
            impl_->plans_[direction].plan, &d_inputGrid, &d_outputGrid, impl_->plans_[direction].info);
}

Gpu3dFft::ImplHipRocfft::ImplHipRocfft(bool                 allocateRealGrid,
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
    Gpu3dFft::Impl::Impl(performOutOfPlaceFFT),
    impl_(std::make_unique<Impl>(allocateRealGrid,
                                 comm,
                                 gridSizesInXForEachRank,
                                 gridSizesInYForEachRank,
                                 nz,
                                 performOutOfPlaceFFT,
                                 context,
                                 pmeStream,
                                 realGridSize,
                                 realGridSizePadded,
                                 complexGridSizePadded,
                                 realGrid,
                                 complexGrid))
{
    allocateComplexGrid(complexGridSizePadded, realGrid, complexGrid, context);
}

Gpu3dFft::ImplHipRocfft::~ImplHipRocfft()
{
    deallocateComplexGrid();
}

} // namespace gmx
