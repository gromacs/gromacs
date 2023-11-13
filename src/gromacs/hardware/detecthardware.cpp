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
#include "gmxpre.h"

#include "detecthardware.h"

#include "config.h"

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/hardware/simd_support.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"

#include "architecture.h"
#include "device_information.h"

#ifdef HAVE_UNISTD_H
#    include <unistd.h> // sysconf()
#endif

gmx_hw_info_t::gmx_hw_info_t(std::unique_ptr<gmx::CpuInfo>          theCpuInfo,
                             std::unique_ptr<gmx::HardwareTopology> theHardwareTopology) :
    cpuInfo(std::move(theCpuInfo)), hardwareTopology(std::move(theHardwareTopology))
{
}

gmx_hw_info_t::~gmx_hw_info_t() = default;

namespace gmx
{

//! Convenience macro to help us avoid ifdefs each time we use sysconf
#if !defined(_SC_NPROCESSORS_ONLN) && defined(_SC_NPROC_ONLN)
#    define _SC_NPROCESSORS_ONLN _SC_NPROC_ONLN
#endif

//! Convenience macro to help us avoid ifdefs each time we use sysconf
#if !defined(_SC_NPROCESSORS_CONF) && defined(_SC_NPROC_CONF)
#    define _SC_NPROCESSORS_CONF _SC_NPROC_CONF
#endif

/*! \brief The result of device detection
 *
 * Note that non-functional device detection still produces
 * a detection result, ie. of no devices. This might not be
 * what the user wanted, so it makes sense to log later when
 * that is possible. */
struct DeviceDetectionResult
{
    //! The device information detected
    std::vector<std::unique_ptr<DeviceInformation>> deviceInfoList_;
    //! Container of possible warnings to issue when that is possible
    std::vector<std::string> deviceDetectionWarnings_;
};

/*! \brief Detect GPUs when that makes sense to attempt.
 *
 * \param[in]  physicalNodeComm  The communicator across this physical node
 * \return The result of the detection, perhaps including diagnostic messages
 *         to issue later.
 *
 * \todo Coordinating the efficient detection of devices across
 * multiple ranks per node should be separated from the lower-level
 * hardware detection. See
 * https://gitlab.com/gromacs/gromacs/-/issues/3650.
 */
static DeviceDetectionResult detectAllDeviceInformation(const PhysicalNodeCommunicator& physicalNodeComm)
{
    DeviceDetectionResult deviceDetectionResult;

    if (!isDeviceDetectionEnabled())
    {
        return deviceDetectionResult;
    }

    std::string errorMessage;

    bool isMainRankOfPhysicalNode = true;
#if GMX_LIB_MPI
    isMainRankOfPhysicalNode = (physicalNodeComm.rank_ == 0);
#else
    // Without an MPI library, this process is trivially the only one
    // on the physical node. This code runs before e.g. thread-MPI
    // ranks are spawned, so detection is race-free by construction.
    // Read-only access is enforced with providing those ranks with a
    // handle to a const object, so usage is also free of races.
    GMX_UNUSED_VALUE(physicalNodeComm);
    isMainRankOfPhysicalNode           = true;
#endif

    /* The SYCL and OpenCL support requires us to run detection on all
     * ranks.
     *
     * With CUDA we don't need to, and prefer to detect on one rank
     * and send the information to the other ranks over MPI. This
     * avoids creating a start-up bottleneck with each MPI rank on a
     * node making the same GPU API calls. */
    constexpr bool allRanksMustDetectGpus = (GMX_GPU_OPENCL != 0 || GMX_GPU_SYCL != 0);
    bool           gpusCanBeDetected      = false;
    if (isMainRankOfPhysicalNode || allRanksMustDetectGpus)
    {
        std::string errorMessage;
        gpusCanBeDetected = isDeviceDetectionFunctional(&errorMessage);
        if (!gpusCanBeDetected)
        {
            deviceDetectionResult.deviceDetectionWarnings_.emplace_back(
                    "Detection of GPUs failed. The API reported:\n" + errorMessage);
        }
    }

    if (gpusCanBeDetected)
    {
        deviceDetectionResult.deviceInfoList_ = findDevices();
        // No need to tell the user anything at this point, they get a
        // hardware report later.
    }

#if GMX_LIB_MPI
    if (!allRanksMustDetectGpus && (physicalNodeComm.size_ > 1))
    {
        // Main rank must serialize the device information list and
        // send it to the other ranks on this node.
        std::vector<char> buffer;
        int               sizeOfBuffer;
        if (isMainRankOfPhysicalNode)
        {
            gmx::InMemorySerializer writer;
            serializeDeviceInformations(deviceDetectionResult.deviceInfoList_, &writer);
            buffer       = writer.finishAndGetBuffer();
            sizeOfBuffer = buffer.size();
        }
        // Ensure all ranks agree on the size of list to be sent
        MPI_Bcast(&sizeOfBuffer, 1, MPI_INT, 0, physicalNodeComm.comm_);
        buffer.resize(sizeOfBuffer);
        if (!buffer.empty())
        {
            // Send the list and deserialize it
            MPI_Bcast(buffer.data(), buffer.size(), MPI_BYTE, 0, physicalNodeComm.comm_);
            if (!isMainRankOfPhysicalNode)
            {
                gmx::InMemoryDeserializer reader(buffer, false);
                deviceDetectionResult.deviceInfoList_ = deserializeDeviceInformations(&reader);
            }
        }
    }
#endif
    return deviceDetectionResult;
}

//! Reduce the locally collected \p hardwareInfo over MPI ranks
static void gmx_collect_hardware_mpi(const gmx::CpuInfo&             cpuInfo,
                                     const PhysicalNodeCommunicator& physicalNodeComm,
                                     gmx_hw_info_t*                  hardwareInfo,
                                     [[maybe_unused]] MPI_Comm       world)
{
    int nCores           = 0;
    int nProcessingUnits = 0;
    for (const auto& p : hardwareInfo->hardwareTopology->machine().packages)
    {
        nCores += p.cores.size();
        for (const auto& c : p.cores)
        {
            nProcessingUnits += c.processingUnits.size();
        }
    }
    int maxThreads = hardwareInfo->hardwareTopology->maxThreads();

    /* Zen1 is assumed for:
     * - family=23 with the below listed models;
     * - Hygon as vendor.
     */
    const bool cpuIsAmdZen1 = gmx::cpuIsAmdZen1(cpuInfo);

    int numCompatibleDevices = getCompatibleDevices(hardwareInfo->deviceInfoList).size();

    // Collect information about GPU-aware MPI support
    const gmx::GpuAwareMpiStatus gpuAwareMpiStatus =
            getMinimalSupportedGpuAwareMpiStatus(hardwareInfo->deviceInfoList);
#if GMX_LIB_MPI
    int gpu_hash;

    /* Create a unique hash of the GPU type(s) in this node */
    gpu_hash = 0;
    /* Here it might be better to only loop over the compatible GPU, but we
     * don't have that information available and it would also require
     * removing the device ID from the device info string.
     */
    for (const auto& deviceInfo : hardwareInfo->deviceInfoList)
    {
        /* Since the device ID is incorporated in the hash, the order of
         * the GPUs affects the hash. Also two identical GPUs won't give
         * a gpu_hash of zero after XORing.
         */
        std::string deviceInfoString = getDeviceInformationString(*deviceInfo);
        gpu_hash ^= gmx_string_fullhash_func(deviceInfoString.c_str(), gmx_string_hash_init);
    }

    constexpr int                      numElementsCounts = 5;
    std::array<int, numElementsCounts> countsReduced;
    {
        std::array<int, numElementsCounts> countsLocal = { { 0 } };
        // Organize to sum values from only one rank within each node,
        // so we get the sum over all nodes.
        bool isMainRankOfPhysicalNode = (physicalNodeComm.rank_ == 0);
        if (isMainRankOfPhysicalNode)
        {
            countsLocal[0] = 1;
            countsLocal[1] = nCores;
            countsLocal[2] = nProcessingUnits;
            countsLocal[3] = maxThreads;
            countsLocal[4] = numCompatibleDevices;
        }

        MPI_Allreduce(countsLocal.data(), countsReduced.data(), countsLocal.size(), MPI_INT, MPI_SUM, world);
    }

    constexpr int                   numElementsMax = 14;
    std::array<int, numElementsMax> maxMinReduced;
    {
        std::array<int, numElementsMax> maxMinLocal;
        /* Store + and - values for all ranks,
         * so we can get max+min with one MPI call.
         */
        maxMinLocal[0]  = nCores;
        maxMinLocal[1]  = nProcessingUnits;
        maxMinLocal[2]  = maxThreads;
        maxMinLocal[3]  = numCompatibleDevices;
        maxMinLocal[4]  = static_cast<int>(gmx::simdSuggested(cpuInfo));
        maxMinLocal[5]  = gpu_hash;
        maxMinLocal[6]  = -maxMinLocal[0];
        maxMinLocal[7]  = -maxMinLocal[1];
        maxMinLocal[8]  = -maxMinLocal[2];
        maxMinLocal[9]  = -maxMinLocal[3];
        maxMinLocal[10] = -maxMinLocal[4];
        maxMinLocal[11] = -maxMinLocal[5];
        maxMinLocal[12] = (cpuIsAmdZen1 ? 1 : 0);
        maxMinLocal[13] =
                -static_cast<int>(gpuAwareMpiStatus); // Enum is ordinal, higher values mean better support

        MPI_Allreduce(maxMinLocal.data(), maxMinReduced.data(), maxMinLocal.size(), MPI_INT, MPI_MAX, world);
    }

    hardwareInfo->nphysicalnode        = countsReduced[0];
    hardwareInfo->ncore_tot            = countsReduced[1];
    hardwareInfo->ncore_min            = -maxMinReduced[6];
    hardwareInfo->ncore_max            = maxMinReduced[0];
    hardwareInfo->nProcessingUnits_tot = countsReduced[2];
    hardwareInfo->nProcessingUnits_min = -maxMinReduced[7];
    hardwareInfo->nProcessingUnits_max = maxMinReduced[1];
    hardwareInfo->maxThreads_tot       = countsReduced[3];
    hardwareInfo->maxThreads_min       = -maxMinReduced[8];
    hardwareInfo->maxThreads_max       = maxMinReduced[2];
    hardwareInfo->ngpu_compatible_tot  = countsReduced[4];
    hardwareInfo->ngpu_compatible_min  = -maxMinReduced[9];
    hardwareInfo->ngpu_compatible_max  = maxMinReduced[3];
    hardwareInfo->simd_suggest_min     = -maxMinReduced[10];
    hardwareInfo->simd_suggest_max     = maxMinReduced[4];
    hardwareInfo->bIdenticalGPUs       = (maxMinReduced[5] == -maxMinReduced[11]);
    hardwareInfo->haveAmdZen1Cpu       = (maxMinReduced[12] > 0);
    hardwareInfo->minGpuAwareMpiStatus = static_cast<gmx::GpuAwareMpiStatus>(-maxMinReduced[13]);

#else
    hardwareInfo->nphysicalnode        = 1;
    hardwareInfo->ncore_tot            = nCores;
    hardwareInfo->ncore_min            = nCores;
    hardwareInfo->ncore_max            = nCores;
    hardwareInfo->nProcessingUnits_tot = nProcessingUnits;
    hardwareInfo->nProcessingUnits_tot = nProcessingUnits;
    hardwareInfo->nProcessingUnits_tot = nProcessingUnits;
    hardwareInfo->maxThreads_tot       = maxThreads;
    hardwareInfo->maxThreads_min       = maxThreads;
    hardwareInfo->maxThreads_max       = maxThreads;
    hardwareInfo->ngpu_compatible_tot  = numCompatibleDevices;
    hardwareInfo->ngpu_compatible_min  = numCompatibleDevices;
    hardwareInfo->ngpu_compatible_max  = numCompatibleDevices;
    hardwareInfo->simd_suggest_min     = static_cast<int>(simdSuggested(cpuInfo));
    hardwareInfo->simd_suggest_max     = static_cast<int>(simdSuggested(cpuInfo));
    hardwareInfo->bIdenticalGPUs       = TRUE;
    hardwareInfo->haveAmdZen1Cpu       = cpuIsAmdZen1;
    hardwareInfo->minGpuAwareMpiStatus = gpuAwareMpiStatus;
    GMX_UNUSED_VALUE(physicalNodeComm);
#endif

    // SIMD instruction sets are actually NOT incremental, even though the
    // code above abuses the enum a bit and assumes that. FMA4 on old AMD
    // CPUs is a particular exception since it is not present on later CPUs,
    // so we at least catch that one explicitly and if there is also a mix
    // nodes we set the lowest supported set to plain 256-bit AVX.
    if (hardwareInfo->simd_suggest_min == static_cast<int>(gmx::SimdType::X86_Avx128Fma)
        && hardwareInfo->simd_suggest_min != hardwareInfo->simd_suggest_max)
    {
        hardwareInfo->simd_suggest_min = static_cast<int>(gmx::SimdType::X86_Avx);
    }
}

std::unique_ptr<gmx_hw_info_t> gmx_detect_hardware(const PhysicalNodeCommunicator& physicalNodeComm,
                                                   MPI_Comm                        libraryCommWorld)
{
    // TODO: We should also do CPU hardware detection only once on each
    // physical node and broadcast it, instead of doing it on every MPI rank.
    auto hardwareInfo = std::make_unique<gmx_hw_info_t>(
            std::make_unique<CpuInfo>(CpuInfo::detect()),
            std::make_unique<HardwareTopology>(HardwareTopology::detect()));

    // Detect GPUs
    // Open a nested scope so no temporary variables can
    // be mis-used later.
    {
        DeviceDetectionResult deviceDetectionResult = detectAllDeviceInformation(physicalNodeComm);
        hardwareInfo->deviceInfoList.swap(deviceDetectionResult.deviceInfoList_);
        std::swap(hardwareInfo->hardwareDetectionWarnings_, deviceDetectionResult.deviceDetectionWarnings_);
    }

    gmx_collect_hardware_mpi(*hardwareInfo->cpuInfo, physicalNodeComm, hardwareInfo.get(), libraryCommWorld);

    return hardwareInfo;
}

void logHardwareDetectionWarnings(const gmx::MDLogger& mdlog, const gmx_hw_info_t& hardwareInformation)
{
    for (const std::string& warningString : hardwareInformation.hardwareDetectionWarnings_)
    {
        GMX_LOG(mdlog.warning).asParagraph().appendText(warningString);
    }
}

} // namespace gmx
