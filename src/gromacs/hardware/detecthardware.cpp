/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016, The GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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
#include "gromacs/simd/support.h"
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
#include "prepare_detection.h"

#ifdef HAVE_UNISTD_H
#    include <unistd.h> // sysconf()
#endif

gmx_hw_info_t::gmx_hw_info_t(std::unique_ptr<gmx::CpuInfo>          cpuInfo,
                             std::unique_ptr<gmx::HardwareTopology> hardwareTopology) :
    cpuInfo(std::move(cpuInfo)),
    hardwareTopology(std::move(hardwareTopology))
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

    bool isMasterRankOfPhysicalNode = true;
#if GMX_LIB_MPI
    isMasterRankOfPhysicalNode = (physicalNodeComm.rank_ == 0);
#else
    // Without an MPI library, this process is trivially the only one
    // on the physical node. This code runs before e.g. thread-MPI
    // ranks are spawned, so detection is race-free by construction.
    // Read-only access is enforced with providing those ranks with a
    // handle to a const object, so usage is also free of races.
    GMX_UNUSED_VALUE(physicalNodeComm);
    isMasterRankOfPhysicalNode        = true;
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
    if (isMasterRankOfPhysicalNode || allRanksMustDetectGpus)
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
        // Master rank must serialize the device information list and
        // send it to the other ranks on this node.
        std::vector<char> buffer;
        int               sizeOfBuffer;
        if (isMasterRankOfPhysicalNode)
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
            if (!isMasterRankOfPhysicalNode)
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
                                     gmx_hw_info_t*                  hardwareInfo)
{
    const int ncore = hardwareInfo->hardwareTopology->numberOfCores();
    /* Zen1 is assumed for:
     * - family=23 with the below listed models;
     * - Hygon as vendor.
     */
    const bool cpuIsAmdZen1 = gmx::cpuIsAmdZen1(cpuInfo);

    int numCompatibleDevices = getCompatibleDevices(hardwareInfo->deviceInfoList).size();
#if GMX_LIB_MPI
    int nhwthread;
    int gpu_hash;

    nhwthread = hardwareInfo->nthreads_hw_avail;
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

    constexpr int                      numElementsCounts = 4;
    std::array<int, numElementsCounts> countsReduced;
    {
        std::array<int, numElementsCounts> countsLocal = { { 0 } };
        // Organize to sum values from only one rank within each node,
        // so we get the sum over all nodes.
        bool isMasterRankOfPhysicalNode = (physicalNodeComm.rank_ == 0);
        if (isMasterRankOfPhysicalNode)
        {
            countsLocal[0] = 1;
            countsLocal[1] = ncore;
            countsLocal[2] = nhwthread;
            countsLocal[3] = numCompatibleDevices;
        }

        MPI_Allreduce(countsLocal.data(), countsReduced.data(), countsLocal.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    constexpr int                   numElementsMax = 11;
    std::array<int, numElementsMax> maxMinReduced;
    {
        std::array<int, numElementsMax> maxMinLocal;
        /* Store + and - values for all ranks,
         * so we can get max+min with one MPI call.
         */
        maxMinLocal[0]  = ncore;
        maxMinLocal[1]  = nhwthread;
        maxMinLocal[2]  = numCompatibleDevices;
        maxMinLocal[3]  = static_cast<int>(gmx::simdSuggested(cpuInfo));
        maxMinLocal[4]  = gpu_hash;
        maxMinLocal[5]  = -maxMinLocal[0];
        maxMinLocal[6]  = -maxMinLocal[1];
        maxMinLocal[7]  = -maxMinLocal[2];
        maxMinLocal[8]  = -maxMinLocal[3];
        maxMinLocal[9]  = -maxMinLocal[4];
        maxMinLocal[10] = (cpuIsAmdZen1 ? 1 : 0);

        MPI_Allreduce(maxMinLocal.data(), maxMinReduced.data(), maxMinLocal.size(), MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }

    hardwareInfo->nphysicalnode       = countsReduced[0];
    hardwareInfo->ncore_tot           = countsReduced[1];
    hardwareInfo->ncore_min           = -maxMinReduced[5];
    hardwareInfo->ncore_max           = maxMinReduced[0];
    hardwareInfo->nhwthread_tot       = countsReduced[2];
    hardwareInfo->nhwthread_min       = -maxMinReduced[6];
    hardwareInfo->nhwthread_max       = maxMinReduced[1];
    hardwareInfo->ngpu_compatible_tot = countsReduced[3];
    hardwareInfo->ngpu_compatible_min = -maxMinReduced[7];
    hardwareInfo->ngpu_compatible_max = maxMinReduced[2];
    hardwareInfo->simd_suggest_min    = -maxMinReduced[8];
    hardwareInfo->simd_suggest_max    = maxMinReduced[3];
    hardwareInfo->bIdenticalGPUs      = (maxMinReduced[4] == -maxMinReduced[9]);
    hardwareInfo->haveAmdZen1Cpu      = (maxMinReduced[10] > 0);
#else
    hardwareInfo->nphysicalnode       = 1;
    hardwareInfo->ncore_tot           = ncore;
    hardwareInfo->ncore_min           = ncore;
    hardwareInfo->ncore_max           = ncore;
    hardwareInfo->nhwthread_tot       = hardwareInfo->nthreads_hw_avail;
    hardwareInfo->nhwthread_min       = hardwareInfo->nthreads_hw_avail;
    hardwareInfo->nhwthread_max       = hardwareInfo->nthreads_hw_avail;
    hardwareInfo->ngpu_compatible_tot = numCompatibleDevices;
    hardwareInfo->ngpu_compatible_min = numCompatibleDevices;
    hardwareInfo->ngpu_compatible_max = numCompatibleDevices;
    hardwareInfo->simd_suggest_min    = static_cast<int>(simdSuggested(cpuInfo));
    hardwareInfo->simd_suggest_max    = static_cast<int>(simdSuggested(cpuInfo));
    hardwareInfo->bIdenticalGPUs      = TRUE;
    hardwareInfo->haveAmdZen1Cpu      = cpuIsAmdZen1;
    GMX_UNUSED_VALUE(physicalNodeComm);
#endif
}

void hardwareTopologyDoubleCheckDetection(const gmx::MDLogger gmx_unused& mdlog,
                                          const gmx::HardwareTopology gmx_unused& hardwareTopology)
{
#if defined HAVE_SYSCONF && defined(_SC_NPROCESSORS_CONF)
    if (hardwareTopology.supportLevel() < gmx::HardwareTopology::SupportLevel::LogicalProcessorCount)
    {
        return;
    }

    int countFromDetection = hardwareTopology.machine().logicalProcessorCount;
    int countConfigured    = sysconf(_SC_NPROCESSORS_CONF);

    /* BIOS, kernel or user actions can take physical processors
     * offline. We already cater for the some of the cases inside the hardwareToplogy
     * by trying to spin up cores just before we detect, but there could be other
     * cases where it is worthwhile to hint that there might be more resources available.
     */
    if (countConfigured >= 0 && countConfigured != countFromDetection)
    {
        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "Note: %d CPUs configured, but only %d were detected to be online.\n",
                        countConfigured,
                        countFromDetection);

        if (c_architecture == Architecture::X86 && countConfigured == 2 * countFromDetection)
        {
            GMX_LOG(mdlog.info)
                    .appendText(
                            "      X86 Hyperthreading is likely disabled; enable it for better "
                            "performance.");
        }
        // For PowerPC (likely Power8) it is possible to set SMT to either 2,4, or 8-way hardware threads.
        // We only warn if it is completely disabled since default performance drops with SMT8.
        if (c_architecture == Architecture::PowerPC && countConfigured == 8 * countFromDetection)
        {
            GMX_LOG(mdlog.info)
                    .appendText(
                            "      PowerPC SMT is likely disabled; enable SMT2/SMT4 for better "
                            "performance.");
        }
    }
#else
    GMX_UNUSED_VALUE(mdlog);
    GMX_UNUSED_VALUE(hardwareTopology);
#endif
}

std::unique_ptr<gmx_hw_info_t> gmx_detect_hardware(const PhysicalNodeCommunicator& physicalNodeComm)
{
    // Ensure all cores have spun up, where applicable.
    hardwareTopologyPrepareDetection();

    // TODO: We should also do CPU hardware detection only once on each
    // physical node and broadcast it, instead of doing it on every MPI rank.
    auto hardwareInfo = std::make_unique<gmx_hw_info_t>(
            std::make_unique<CpuInfo>(CpuInfo::detect()),
            std::make_unique<HardwareTopology>(HardwareTopology::detect()));

    // TODO: Get rid of this altogether.
    hardwareInfo->nthreads_hw_avail = hardwareInfo->hardwareTopology->machine().logicalProcessorCount;

    // Detect GPUs
    // Open a nested scope so no temporary variables can
    // be mis-used later.
    {
        DeviceDetectionResult deviceDetectionResult = detectAllDeviceInformation(physicalNodeComm);
        hardwareInfo->deviceInfoList.swap(deviceDetectionResult.deviceInfoList_);
        std::swap(hardwareInfo->hardwareDetectionWarnings_, deviceDetectionResult.deviceDetectionWarnings_);
    }

    gmx_collect_hardware_mpi(*hardwareInfo->cpuInfo, physicalNodeComm, hardwareInfo.get());

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
