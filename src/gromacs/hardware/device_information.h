/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
/*! \libinternal \file
 *  \brief Declares the GPU information structure and its helpers
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef GMX_HARDWARE_DEVICE_INFORMATION_H
#define GMX_HARDWARE_DEVICE_INFORMATION_H

#include "config.h"

#include <array>
#include <optional>
#include <type_traits>

#if GMX_GPU_CUDA
#    include <cuda_runtime.h>
#endif

#if GMX_GPU_HIP
#    include <hip/hip_runtime.h>
#endif

#if GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/gmxopencl.h"
#endif

#if GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gmxsycl.h"
#endif

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fixedcapacityvector.h"
#include "gromacs/utility/mpiinfo.h"

//! Constant used to help minimize preprocessed code
static constexpr bool c_binarySupportsGpus = (GMX_GPU != 0);

//! Possible results of the GPU detection/check.
enum class DeviceStatus : int
{
    //! The device is compatible
    Compatible,
    //! Device does not exist
    Nonexistent,
    //! Device is not compatible
    Incompatible,
    //! OpenCL device has incompatible cluster size for non-bonded kernels.
    IncompatibleClusterSize,
    //! There are known issues with OpenCL on NVIDIA Volta and newer.
    IncompatibleNvidiaVolta,
    /*! \brief The device originates from non-recommended SYCL backend.
     * The device might work by itself, but to simplify device allocation, it is marked as incompatible.
     * */
    NotPreferredBackend,
    /*! \brief An error occurred during the functionality checks.
     * That indicates malfunctioning of the device, driver, or incompatible driver/runtime.
     */
    NonFunctional,
    /*! \brief CUDA devices are busy or unavailable.
     * typically due to use of \p cudaComputeModeExclusive, \p cudaComputeModeProhibited modes.
     */
    Unavailable,
    /*! \brief The device is outside the set of compilation targets.
     * See \c GMX_CUDA_TARGET_SM and \c GMX_CUDA_TARGET_COMPUTE CMake variables.
     */
    DeviceNotTargeted,
    //! \brief AMD RDNA devices (gfx10xx, gfx11xx) with 32-wide execution are not supported with OpenCL,
    IncompatibleOclAmdRdna,
    //! \brief RDNA not targeted (SYCL)
    IncompatibleAmdRdnaNotTargeted,
    //! Enumeration size
    Count
};

/*! \brief Names of the GPU detection/check results
 *
 * Check-source wants to warn about the use of a symbol name that would
 * require an inclusion of config.h. However the use is in a comment, so that
 * is a false warning. So C-style string concatenation is used to fool the
 * naive parser in check-source. That needs a clang-format suppression
 * in order to look reasonable. Also clang-tidy wants to suggest that a comma is
 * missing, so that is suppressed.
 */
static const gmx::EnumerationArray<DeviceStatus, const char*> c_deviceStateString = {
    "compatible",
    "nonexistent",
    "incompatible",
    // clang-format off
    // NOLINTNEXTLINE(bugprone-suspicious-missing-comma)
    "incompatible (please recompile with correct GMX" "_GPU_NB_CLUSTER_SIZE of 4)",
    // clang-format on
    "incompatible (please use CUDA build for NVIDIA Volta GPUs or newer)",
    "not recommended (please use SYCL_DEVICE_FILTER to limit visibility to a single backend)",
    "non-functional",
    "unavailable",
    "not in set of targeted devices",
    "incompatible (AMD RDNA devices are not supported)", // Issue #4521
    // clang-format off
    // NOLINTNEXTLINE(bugprone-suspicious-missing-comma)
    "incompatible (please recompile with GMX" "_ACPP_ENABLE_AMD_RDNA_SUPPORT)"
    // clang-format on
};

//! Device vendors
enum class DeviceVendor : int
{
    //! No data
    Unknown = 0,
    //! NVIDIA
    Nvidia = 1,
    //! Advanced Micro Devices
    Amd = 2,
    //! Intel
    Intel = 3,
    //! Apple
    Apple = 4,
    //! Enumeration size
    Count = 5
};


/*! \libinternal \brief Platform-dependent device information.
 *
 * The device information is queried and set at detection and contains
 * both information about the device/hardware returned by the runtime as well
 * as additional data like support status.
 */
struct DeviceInformation
{
    //! Device status.
    DeviceStatus status;
    //! ID of the device.
    int id;
    //! Device vendor.
    DeviceVendor deviceVendor;
    /*! \brief Warp/sub-group sizes supported by the device.
     *
     * \ref DeviceInformation must be serializable in CUDA, so we cannot use \c std::vector here.
     * Arbitrarily limiting to 10.
     */
    gmx::FixedCapacityVector<int, 10> supportedSubGroupSizes;

    gmx::GpuAwareMpiStatus gpuAwareMpiStatus;
#if GMX_GPU_CUDA
    //! CUDA device properties.
    cudaDeviceProp prop;
#elif GMX_GPU_HIP
    //! HIP device properties.
    hipDeviceProp_t prop;
#elif GMX_GPU_OPENCL
    cl_platform_id oclPlatformId;       //!< OpenCL Platform ID.
    cl_device_id   oclDeviceId;         //!< OpenCL Device ID.
    char           device_name[256];    //!< Device name.
    char           device_version[256]; //!< Device version.
    char           vendorName[256];     //!< Device vendor name.
    int            compute_units;       //!< Number of compute units.
    int            adress_bits;         //!< Number of address bits the device is capable of.
    size_t         maxWorkItemSizes[3]; //!< Workgroup size limits (CL_DEVICE_MAX_WORK_ITEM_SIZES).
    size_t         maxWorkGroupSize;    //!< Workgroup total size limit (CL_DEVICE_MAX_WORK_GROUP_SIZE).
#elif GMX_GPU_SYCL
    sycl::device syclDevice;
    //! CUDA CC major for NVIDIA devices, generation code for AMD (gfx90a -> 9), not set for Intel (yet)
    std::optional<int> hardwareVersionMajor;
    //! CUDA CC minor for NVIDIA devices, major architecture(?) code for AMD (gfx90a -> 0), not set for Intel (yet)
    std::optional<int> hardwareVersionMinor;
    //! CUDA CC minor for NVIDIA devices, device code for AMD (gfx90a -> a -> 10), not set for Intel (yet)
    std::optional<int> hardwareVersionPatch;
#endif
};

//! Whether \ref DeviceInformation can be serialized for sending via MPI.
static constexpr bool c_canSerializeDeviceInformation = std::is_trivially_copyable_v<DeviceInformation>;

#endif // GMX_HARDWARE_DEVICE_INFORMATION_H
