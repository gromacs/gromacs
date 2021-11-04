/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 *  \brief Implements GPU 3D FFT routines for hipSYCL via rocFFT.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * For hipSYCL, in order to call FFT APIs from the respective vendors
 * using the same DeviceStream as other operations, a vendor extension
 * called "custom operations" is used (see hipSYCL
 * doc/enqueue-custom-operation.md). That effectively enqueues an
 * asynchronous host-side lambda into the same queue. The body of the
 * lambda unpacks the runtime data structures to get the native
 * handles and calls the native FFT APIs.
 *
 * For a 3D FFT, rocFFT requires a working buffer which it allocates
 * itself if not provided. This might be slow enough to be worth
 * optimizing. This working buffer could be provided in advance by
 * calling rocfft_plan_get_work_buffer_size, allocating a buffer that
 * persists suitably, and then using
 * rocfft_execution_info_set_work_buffer in a custom operation.
 *
 * hipSYCL queues operate at a higher level of abstraction than hip
 * streams, with the runtime distributing work to the latter to
 * balance load. It is possible to set the HIP stream in
 * rocfft_execution_info, but then there is no guarantee that a
 * subsequent queue item will run using the same stream. So we
 * currently do not attempt to set the stream.
 *
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_sycl_rocfft.h"

#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"

#include <vector>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#ifndef __HIPSYCL__
#    error This file can only be compiled with hipSYCL enabled
#endif

#if !defined HIPSYCL_PLATFORM_ROCM
#    error Only ROCM platform is supported for 3D FFT with hipSYCL
#endif

#include "rocfft.h"

namespace gmx
{

namespace
{

//! Model the kinds of 3D FFT implemented
enum class FftDirection : int
{
    RealToComplex,
    ComplexToReal,
    Count,
};

//! Strings that match enum rocfft_status_e in rocfft.h
const std::array<const char*, rocfft_status_invalid_work_buffer + 1> c_rocfftErrorStrings = {
    "success",
    "failure",
    "invalid argument value",
    "invalid dimensions",
    "invalid array type",
    "invalid strides",
    "invalid distance",
    "invalid offset",
    "invalid work buffer"
};

//! Helper for consistent error handling
void handleFftError(rocfft_status result, const std::string& msg)
{
    if (result != rocfft_status_success)
    {
        if (result <= rocfft_status_invalid_work_buffer)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString(
                    "%s: (error code %d - %s)\n", msg.c_str(), result, c_rocfftErrorStrings[result])));
        }
        else
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("%s: (error code %d)\n", msg.c_str(), result)));
        }
    }
}

//! Helper for consistent error handling
void handleFftError(rocfft_status result, const std::string& direction, const std::string& msg)
{
    if (result != rocfft_status_success)
    {
        handleFftError(result, msg + " doing " + direction);
    }
}

//! Provides RAII-style initialization of rocFFT library
class RocfftInitializer
{
public:
    RocfftInitializer()
    {
        rocfft_status result;
        result = rocfft_setup();
        handleFftError(result, "rocfft_setup failure");
    }
    ~RocfftInitializer()
    {
        // No need to handle any errors in a destructor, and
        // anyway one cannot throw.
        rocfft_cleanup();
    }
};

//! All the persistent data for planning an executing a 3D FFT
struct RocfftPlan
{
    //! Describes details of the data layout
    rocfft_plan_description description = nullptr;
    //! High level information about the plan
    rocfft_plan plan = nullptr;
    //! Destructor
    ~RocfftPlan()
    {
        // No need to handle any errors in a destructor,
        // and anyway one cannot throw.
        if (plan)
        {
            rocfft_plan_destroy(plan);
        }
        if (description)
        {
            rocfft_plan_description_destroy(description);
        }
    }
};

//! Helper struct to reduce repetitive code setting up a 3D FFT plan
struct PlanSetupData
{
    //! Format of the input array (real or hermitian)
    rocfft_array_type arrayType;
    //! Strides through the input array for the three dimensions
    std::array<size_t, DIM> strides;
    //! Total size of the input array (including padding)
    size_t totalSize;
};

//! Compute the stride through the real 1D array
std::array<size_t, DIM> makeRealStrides(ivec realGridSizePadded)
{
    return { 1, size_t(realGridSizePadded[ZZ]), size_t(realGridSizePadded[ZZ] * realGridSizePadded[YY]) };
};

//! Compute the stride through the complex 1D array
std::array<size_t, DIM> makeComplexStrides(ivec complexGridSizePadded)
{
    return { 1,
             size_t(complexGridSizePadded[XX]),
             size_t(complexGridSizePadded[XX] * complexGridSizePadded[YY]) };
}

//! Compute total grid size
size_t computeTotalSize(ivec gridSize)
{
    return size_t(gridSize[XX] * gridSize[YY] * gridSize[ZZ]);
}

/*! \brief Prepare plans for the forward and reverse transformation.
 *
 * Because these require device-side allocations, some of them must be
 * done in a SYCL queue. */
RocfftPlan makePlan(const std::string&     descriptiveString,
                    rocfft_transform_type  transformType,
                    const PlanSetupData&   inputPlanSetupData,
                    const PlanSetupData&   outputPlanSetupData,
                    ArrayRef<const size_t> rocfftRealGridSize,
                    const DeviceStream&    pmeStream)
{
    rocfft_plan_description description = nullptr;
    rocfft_status           result;
    result = rocfft_plan_description_create(&description);
    handleFftError(result, descriptiveString, "rocfft_plan_description_create failure");
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
    handleFftError(result, descriptiveString, "rocfft_plan_description_set_data_layout failure");

    // The plan creation depends on the identity of the GPU device, so
    // we make sure it is made in the same queue where it will be
    // used. The stream for execution can be set at the same time.

    // First set up device buffers to receive the rocfft status values
    rocfft_plan                        plan = nullptr;
    cl::sycl::buffer<rocfft_status, 1> resultPlanCreate(1);

    // Submit the planning to the queue. This is necessary so that we
    // can ensure that the allocations in the planning go to the right
    // context.
    {
        auto queue = pmeStream.stream();
        // Make a buffer that is a view of the existing memory for a
        // plan.
        cl::sycl::buffer<rocfft_plan, 1> planView =
                cl::sycl::make_async_writeback_view(&plan, cl::sycl::range(1), queue);
        queue.submit([&](cl::sycl::handler& cgh) {
            // Make the necessary accessors
            auto a_plan = planView.get_access(cgh, cl::sycl::write_only, cl::sycl::no_init);
            auto a_resultPlanCreate =
                    resultPlanCreate.get_access(cgh, cl::sycl::write_only, cl::sycl::no_init);
            cgh.hipSYCL_enqueue_custom_operation([=](cl::sycl::interop_handle& /*h*/) {
                const int numBatches = 1;
                // Unlike some other FFT APIs, in rocFFT the
                // dimension of an FFT is the (vectorial) size
                // of the problem (ie. rocfftRealGridSize), not
                // the size of the input vector (which varies
                // according to whether the input format is real
                // or hermitian respectively for forward or
                // reverse transforms).
                a_resultPlanCreate[0] = rocfft_plan_create(&a_plan[0],
                                                           rocfft_placement_notinplace,
                                                           transformType,
                                                           rocfft_precision_single,
                                                           rocfftRealGridSize.size(),
                                                           rocfftRealGridSize.data(),
                                                           numBatches,
                                                           description);
            });
        });
    }
    // Check for errors that happened while running the hipSYCL custom
    // operation.
    handleFftError(
            resultPlanCreate.get_host_access()[0], descriptiveString, "rocfft_plan_create failure");

    return RocfftPlan{ description, plan };
}

} // namespace

//! Impl class
class Gpu3dFft::ImplSyclRocfft::Impl
{
public:
    //! \copydoc Gpu3dFft::Impl::Impl
    Impl(bool                 allocateGrids,
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
    cl::sycl::buffer<float, 1> realGrid_;
    //! Handle to the complex grid buffer
    cl::sycl::buffer<float, 1> complexGrid_;
    /*! \brief Copy of PME stream
     *
     * This copy is guaranteed by the SYCL standard to work as if
     * it was the original. */
    cl::sycl::queue queue_;
};

Gpu3dFft::ImplSyclRocfft::Impl::Impl(bool allocateGrids,
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
                                     DeviceBuffer<float>* complexGrid) :
    plans_{
        makePlan("real-to-complex",
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
                                      size_t(realGridSize[XX]) },
                 pmeStream),
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
                                      size_t(realGridSize[XX]) },
                 pmeStream),
    },
    realGrid_(*realGrid->buffer_.get()),
    complexGrid_(*complexGrid->buffer_.get()),
    queue_(pmeStream.stream())
{
    GMX_RELEASE_ASSERT(performOutOfPlaceFFT, "Only out-of-place FFT is implemented in hipSYCL");
    GMX_RELEASE_ASSERT(allocateGrids == false, "Grids need to be pre-allocated");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() == 1 && gridSizesInYForEachRank.size() == 1,
                       "FFT decomposition not implemented with SYCL backend");
}

void Gpu3dFft::ImplSyclRocfft::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
    GMX_RELEASE_ASSERT((dir == GMX_FFT_REAL_TO_COMPLEX) || (dir == GMX_FFT_COMPLEX_TO_REAL),
                       "Only real-to-complex and complex-to-real FFTs are implemented in hipSYCL");
    FftDirection               direction;
    cl::sycl::buffer<float, 1>*inputGrid = nullptr, *outputGrid = nullptr;
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
    // Enqueue the 3D FFT work
    impl_->queue_.submit([&](cl::sycl::handler& cgh) {
        auto inputGridAccessor = inputGrid->get_access(cgh, cl::sycl::read_only);
        auto outputGridAccessor = outputGrid->get_access(cgh, cl::sycl::write_only, cl::sycl::no_init);
        // Use a hipSYCL custom operation to access the native buffers
        // needed to call rocFFT
        cgh.hipSYCL_enqueue_custom_operation([=](cl::sycl::interop_handle& h) {
            void* d_inputGrid  = h.get_native_mem<cl::sycl::backend::hip>(inputGridAccessor);
            void* d_outputGrid = h.get_native_mem<cl::sycl::backend::hip>(outputGridAccessor);
            // Don't check results generated asynchronously,
            // because we don't know what to do with them
            rocfft_execute(impl_->plans_[direction].plan, &d_inputGrid, &d_outputGrid, nullptr);
        });
    });
}

Gpu3dFft::ImplSyclRocfft::ImplSyclRocfft(bool                 allocateGrids,
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
    impl_(std::make_unique<Impl>(allocateGrids,
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
}

Gpu3dFft::ImplSyclRocfft::~ImplSyclRocfft() = default;

} // namespace gmx
