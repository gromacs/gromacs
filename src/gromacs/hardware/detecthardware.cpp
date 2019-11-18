/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#include <chrono>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "gromacs/compat/pointers.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/cpuinfo.h"
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
#include "gromacs/utility/logger.h"
#include "gromacs/utility/mutex.h"
#include "gromacs/utility/physicalnodecommunicator.h"

#include "architecture.h"

#ifdef HAVE_UNISTD_H
#    include <unistd.h> // sysconf()
#endif

gmx_hw_info_t::gmx_hw_info_t(std::unique_ptr<gmx::CpuInfo>          cpuInfo,
                             std::unique_ptr<gmx::HardwareTopology> hardwareTopology) :
    cpuInfo(std::move(cpuInfo)),
    hardwareTopology(std::move(hardwareTopology))
{
}

gmx_hw_info_t::~gmx_hw_info_t()
{
    free_gpu_info(&gpu_info);
}

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

/*! \brief Information about the hardware of all nodes (common to all threads in this process).
 *
 * This information is constructed only when required, but thereafter
 * its lifetime is that of the whole process, potentially across
 * multiple successive simulation parts. It's wise to ensure that only
 * one thread can create the information, but thereafter they can all
 * read it without e.g. needing a std::shared_ptr to ensure its
 * lifetime exceeds that of the thread. */
static std::unique_ptr<gmx_hw_info_t> g_hardwareInfo;
//! A mutex to protect the hwinfo structure
static Mutex g_hardwareInfoMutex;

//! Detect GPUs, if that makes sense to attempt.
static void gmx_detect_gpus(const gmx::MDLogger&             mdlog,
                            const PhysicalNodeCommunicator&  physicalNodeComm,
                            compat::not_null<gmx_hw_info_t*> hardwareInfo)
{
    hardwareInfo->gpu_info.bDetectGPUs = canPerformGpuDetection();

    if (!hardwareInfo->gpu_info.bDetectGPUs)
    {
        return;
    }

    bool isMasterRankOfPhysicalNode = true;
#if GMX_LIB_MPI
    isMasterRankOfPhysicalNode = (physicalNodeComm.rank_ == 0);
#else
    // We choose to run the detection only once with thread-MPI and
    // use a mutex to enforce it.
    GMX_UNUSED_VALUE(physicalNodeComm);
    isMasterRankOfPhysicalNode = true;
#endif

    /* The OpenCL support requires us to run detection on all ranks.
     * With CUDA we don't need to, and prefer to detect on one rank
     * and send the information to the other ranks over MPI. */
    bool allRanksMustDetectGpus = (GMX_GPU == GMX_GPU_OPENCL);
    bool gpusCanBeDetected      = false;
    if (isMasterRankOfPhysicalNode || allRanksMustDetectGpus)
    {
        std::string errorMessage;
        gpusCanBeDetected = isGpuDetectionFunctional(&errorMessage);
        if (!gpusCanBeDetected)
        {
            GMX_LOG(mdlog.info)
                    .asParagraph()
                    .appendTextFormatted(
                            "NOTE: Detection of GPUs failed. The API reported:\n"
                            "      %s\n"
                            "      GROMACS cannot run tasks on a GPU.",
                            errorMessage.c_str());
        }
    }

    if (gpusCanBeDetected)
    {
        findGpus(&hardwareInfo->gpu_info);
        // No need to tell the user anything at this point, they get a
        // hardware report later.
    }

#if GMX_LIB_MPI
    if (!allRanksMustDetectGpus)
    {
        /* Broadcast the GPU info to the other ranks within this node */
        MPI_Bcast(&hardwareInfo->gpu_info.n_dev, 1, MPI_INT, 0, physicalNodeComm.comm_);

        if (hardwareInfo->gpu_info.n_dev > 0)
        {
            int dev_size;

            dev_size = hardwareInfo->gpu_info.n_dev * sizeof_gpu_dev_info();

            if (!isMasterRankOfPhysicalNode)
            {
                hardwareInfo->gpu_info.gpu_dev = (struct gmx_device_info_t*)malloc(dev_size);
            }
            MPI_Bcast(hardwareInfo->gpu_info.gpu_dev, dev_size, MPI_BYTE, 0, physicalNodeComm.comm_);
            MPI_Bcast(&hardwareInfo->gpu_info.n_dev_compatible, 1, MPI_INT, 0, physicalNodeComm.comm_);
        }
    }
#endif
}

//! Reduce the locally collected \p hardwareInfo over MPI ranks
static void gmx_collect_hardware_mpi(const gmx::CpuInfo&              cpuInfo,
                                     const PhysicalNodeCommunicator&  physicalNodeComm,
                                     compat::not_null<gmx_hw_info_t*> hardwareInfo)
{
    const int ncore = hardwareInfo->hardwareTopology->numberOfCores();
    /* Zen1 is assumed for:
     * - family=23 with the below listed models;
     * - Hygon as vendor.
     */
    const bool cpuIsAmdZen1 = ((cpuInfo.vendor() == CpuInfo::Vendor::Amd && cpuInfo.family() == 23
                                && (cpuInfo.model() == 1 || cpuInfo.model() == 17
                                    || cpuInfo.model() == 8 || cpuInfo.model() == 24))
                               || cpuInfo.vendor() == CpuInfo::Vendor::Hygon);
#if GMX_LIB_MPI
    int nhwthread, ngpu, i;
    int gpu_hash;

    nhwthread = hardwareInfo->nthreads_hw_avail;
    ngpu      = hardwareInfo->gpu_info.n_dev_compatible;
    /* Create a unique hash of the GPU type(s) in this node */
    gpu_hash = 0;
    /* Here it might be better to only loop over the compatible GPU, but we
     * don't have that information available and it would also require
     * removing the device ID from the device info string.
     */
    for (i = 0; i < hardwareInfo->gpu_info.n_dev; i++)
    {
        char stmp[STRLEN];

        /* Since the device ID is incorporated in the hash, the order of
         * the GPUs affects the hash. Also two identical GPUs won't give
         * a gpu_hash of zero after XORing.
         */
        get_gpu_device_info_string(stmp, hardwareInfo->gpu_info, i);
        gpu_hash ^= gmx_string_fullhash_func(stmp, gmx_string_hash_init);
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
            countsLocal[3] = ngpu;
        }

        MPI_Allreduce(countsLocal.data(), countsReduced.data(), countsLocal.size(), MPI_INT,
                      MPI_SUM, MPI_COMM_WORLD);
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
        maxMinLocal[2]  = ngpu;
        maxMinLocal[3]  = static_cast<int>(gmx::simdSuggested(cpuInfo));
        maxMinLocal[4]  = gpu_hash;
        maxMinLocal[5]  = -maxMinLocal[0];
        maxMinLocal[6]  = -maxMinLocal[1];
        maxMinLocal[7]  = -maxMinLocal[2];
        maxMinLocal[8]  = -maxMinLocal[3];
        maxMinLocal[9]  = -maxMinLocal[4];
        maxMinLocal[10] = (cpuIsAmdZen1 ? 1 : 0);

        MPI_Allreduce(maxMinLocal.data(), maxMinReduced.data(), maxMinLocal.size(), MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD);
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
    /* All ranks use the same pointer, protected by a mutex in the caller */
    hardwareInfo->nphysicalnode       = 1;
    hardwareInfo->ncore_tot           = ncore;
    hardwareInfo->ncore_min           = ncore;
    hardwareInfo->ncore_max           = ncore;
    hardwareInfo->nhwthread_tot       = hardwareInfo->nthreads_hw_avail;
    hardwareInfo->nhwthread_min       = hardwareInfo->nthreads_hw_avail;
    hardwareInfo->nhwthread_max       = hardwareInfo->nthreads_hw_avail;
    hardwareInfo->ngpu_compatible_tot = hardwareInfo->gpu_info.n_dev_compatible;
    hardwareInfo->ngpu_compatible_min = hardwareInfo->gpu_info.n_dev_compatible;
    hardwareInfo->ngpu_compatible_max = hardwareInfo->gpu_info.n_dev_compatible;
    hardwareInfo->simd_suggest_min    = static_cast<int>(simdSuggested(cpuInfo));
    hardwareInfo->simd_suggest_max    = static_cast<int>(simdSuggested(cpuInfo));
    hardwareInfo->bIdenticalGPUs      = TRUE;
    hardwareInfo->haveAmdZen1Cpu      = cpuIsAmdZen1;
    GMX_UNUSED_VALUE(physicalNodeComm);
#endif
}

/*! \brief Utility that does dummy computing for max 2 seconds to spin up cores
 *
 *  This routine will check the number of cores configured and online
 *  (using sysconf), and the spins doing dummy compute operations for up to
 *  2 seconds, or until all cores have come online. This can be used prior to
 *  hardware detection for platforms that take unused processors offline.
 *
 *  This routine will not throw exceptions.
 */
static void spinUpCore() noexcept
{
#if defined(HAVE_SYSCONF) && defined(_SC_NPROCESSORS_CONF) && defined(_SC_NPROCESSORS_ONLN)
    float dummy           = 0.1;
    int   countConfigured = sysconf(_SC_NPROCESSORS_CONF);    // noexcept
    auto  start           = std::chrono::steady_clock::now(); // noexcept

    while (sysconf(_SC_NPROCESSORS_ONLN) < countConfigured
           && std::chrono::steady_clock::now() - start < std::chrono::seconds(2))
    {
        for (int i = 1; i < 10000; i++)
        {
            dummy /= i;
        }
    }

    if (dummy < 0)
    {
        printf("This cannot happen, but prevents loop from being optimized away.");
    }
#endif
}

/*! \brief Prepare the system before hardware topology detection
 *
 * This routine should perform any actions we want to put the system in a state
 * where we want it to be before detecting the hardware topology. For most
 * processors there is nothing to do, but some architectures (in particular ARM)
 * have support for taking configured cores offline, which will make them disappear
 * from the online processor count.
 *
 * This routine checks if there is a mismatch between the number of cores
 * configured and online, and in that case we issue a small workload that
 * attempts to wake sleeping cores before doing the actual detection.
 *
 * This type of mismatch can also occur for x86 or PowerPC on Linux, if SMT has only
 * been disabled in the kernel (rather than bios). Since those cores will never
 * come online automatically, we currently skip this test for x86 & PowerPC to
 * avoid wasting 2 seconds. We also skip the test if there is no thread support.
 *
 * \note Cores will sleep relatively quickly again, so it's important to issue
 *       the real detection code directly after this routine.
 */
static void hardwareTopologyPrepareDetection()
{
#if defined(HAVE_SYSCONF) && defined(_SC_NPROCESSORS_CONF) \
        && (defined(THREAD_PTHREADS) || defined(THREAD_WINDOWS))

    // Modify this conditional when/if x86 or PowerPC starts to sleep some cores
    if (c_architecture != Architecture::X86 && c_architecture != Architecture::PowerPC)
    {
        int                      countConfigured = sysconf(_SC_NPROCESSORS_CONF);
        std::vector<std::thread> workThreads(countConfigured);

        for (auto& t : workThreads)
        {
            t = std::thread(spinUpCore);
        }

        for (auto& t : workThreads)
        {
            t.join();
        }
    }
#endif
}

/*! \brief Sanity check hardware topology and print some notes to log
 *
 *  \param mdlog            Logger.
 *  \param hardwareTopology Reference to hardwareTopology object.
 */
static void hardwareTopologyDoubleCheckDetection(const gmx::MDLogger gmx_unused& mdlog,
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
                        countConfigured, countFromDetection);

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
#endif
}

gmx_hw_info_t* gmx_detect_hardware(const gmx::MDLogger& mdlog, const PhysicalNodeCommunicator& physicalNodeComm)
{
    // By construction, only one thread ever runs hardware detection,
    // but we may as well prevent issues arising if that would change.
    // Taking the lock early ensures that exactly one thread can
    // attempt to construct g_hardwareInfo.
    lock_guard<Mutex> lock(g_hardwareInfoMutex);

    // If we already have the information, just return a handle to it.
    if (g_hardwareInfo != nullptr)
    {
        return g_hardwareInfo.get();
    }

    // Make the new hardwareInfo in a temporary.
    hardwareTopologyPrepareDetection();

    // TODO: We should also do CPU hardware detection only once on each
    // physical node and broadcast it, instead of doing it on every MPI rank.
    auto hardwareInfo = std::make_unique<gmx_hw_info_t>(
            std::make_unique<CpuInfo>(CpuInfo::detect()),
            std::make_unique<HardwareTopology>(HardwareTopology::detect()));

    // If we detected the topology on this system, double-check that it makes sense
    if (hardwareInfo->hardwareTopology->isThisSystem())
    {
        hardwareTopologyDoubleCheckDetection(mdlog, *hardwareInfo->hardwareTopology);
    }

    // TODO: Get rid of this altogether.
    hardwareInfo->nthreads_hw_avail = hardwareInfo->hardwareTopology->machine().logicalProcessorCount;

    // Detect GPUs
    hardwareInfo->gpu_info.n_dev            = 0;
    hardwareInfo->gpu_info.n_dev_compatible = 0;
    hardwareInfo->gpu_info.gpu_dev          = nullptr;

    gmx_detect_gpus(mdlog, physicalNodeComm, compat::make_not_null(hardwareInfo));
    gmx_collect_hardware_mpi(*hardwareInfo->cpuInfo, physicalNodeComm, compat::make_not_null(hardwareInfo));

    // Now that the temporary is fully constructed, swap it to become
    // the real thing.
    g_hardwareInfo.swap(hardwareInfo);

    return g_hardwareInfo.get();
}

bool compatibleGpusFound(const gmx_gpu_info_t& gpu_info)
{
    return gpu_info.n_dev_compatible > 0;
}

} // namespace gmx
