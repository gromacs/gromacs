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

#include "printhardware.h"

#include "config.h"

#include <cstdlib>

#include <filesystem>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/hardware/identifyavx512fmaunits.h"
#include "gromacs/hardware/simd_support.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#include "architecture.h"

struct DeviceInformation;

//! Constant used to help minimize preprocessed code. Currently disabled also for HIP builds until device detection is implemented.
static constexpr bool bGPUBinary = (GMX_GPU != 0) && !GMX_GPU_HIP;

/*! \internal \brief
 * Returns the GPU information text, one GPU per line.
 */
static std::string sprint_gpus(const std::vector<std::unique_ptr<DeviceInformation>>& deviceInfoList)
{
    std::vector<std::string> gpuStrings(0);
    for (const auto& deviceInfo : deviceInfoList)
    {
        gpuStrings.emplace_back("    " + getDeviceInformationString(*deviceInfo));
    }
    return gmx::joinStrings(gpuStrings, "\n");
}

/* Give a suitable fatal error or warning if the build configuration
   and runtime CPU do not match. */
static void check_use_of_rdtscp_on_this_cpu(const gmx::MDLogger& mdlog, const gmx::CpuInfo& cpuInfo)
{
    bool binaryUsesRdtscp = GMX_USE_RDTSCP;

    const char* programName = gmx::getProgramContext().displayName();

    if (cpuInfo.supportLevel() < gmx::CpuInfo::SupportLevel::Features)
    {
        if (binaryUsesRdtscp)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "The %s executable was compiled to use the rdtscp CPU instruction. "
                            "We cannot detect the features of your current CPU, but will proceed "
                            "anyway. "
                            "If you get a crash, rebuild GROMACS with the GMX_USE_RDTSCP=OFF CMake "
                            "option.",
                            programName);
        }
    }
    else
    {
        bool cpuHasRdtscp = cpuInfo.feature(gmx::CpuInfo::Feature::X86_Rdtscp);

        if (!cpuHasRdtscp && binaryUsesRdtscp)
        {
            gmx_fatal(FARGS,
                      "The %s executable was compiled to use the rdtscp CPU instruction. "
                      "However, this is not supported by the current hardware and continuing would "
                      "lead to a crash. "
                      "Please rebuild GROMACS with the GMX_USE_RDTSCP=OFF CMake option.",
                      programName);
        }

        if (cpuHasRdtscp && !binaryUsesRdtscp)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "The current CPU can measure timings more accurately than the code in\n"
                            "%s was configured to use. This might affect your simulation\n"
                            "speed as accurate timings are needed for load-balancing.\n"
                            "Please consider rebuilding %s with the GMX_USE_RDTSCP=ON CMake "
                            "option.",
                            programName,
                            programName);
        }
    }
}

static std::string detected_hardware_string(const gmx_hw_info_t* hwinfo, bool bFullCpuInfo)
{
    std::string s;

    const gmx::CpuInfo&          cpuInfo = *hwinfo->cpuInfo;
    const gmx::HardwareTopology& hwTop   = *hwinfo->hardwareTopology;

    s = gmx::formatString("\n");
    s += gmx::formatString("Running on %d node%s with total",
                           hwinfo->nphysicalnode,
                           hwinfo->nphysicalnode == 1 ? "" : "s");
    if (hwinfo->ncore_tot > 0)
    {
        s += gmx::formatString(" %d cores,", hwinfo->ncore_tot);
    }
    s += gmx::formatString(" %d processing units", hwinfo->nProcessingUnits_tot);
    if (canPerformDeviceDetection(nullptr))
    {
        s += gmx::formatString(", %d compatible GPU%s",
                               hwinfo->ngpu_compatible_tot,
                               hwinfo->ngpu_compatible_tot == 1 ? "" : "s");
    }
    else if (bGPUBinary)
    {
        if (isDeviceDetectionEnabled())
        {
            s += gmx::formatString(" (GPU detection failed)");
        }
        else
        {
            s += gmx::formatString(" (GPU detection deactivated)");
        }
    }
    s += gmx::formatString("\n");

    if (hwinfo->nphysicalnode > 1)
    {
        /* Print per node hardware feature counts */
        if (hwinfo->ncore_max > 0)
        {
            s += gmx::formatString("  Cores per node:           %2d", hwinfo->ncore_min);
            if (hwinfo->ncore_max > hwinfo->ncore_min)
            {
                s += gmx::formatString(" - %2d", hwinfo->ncore_max);
            }
            s += gmx::formatString("\n");
        }
        s += gmx::formatString("  Logical processing units per node:   %2d", hwinfo->nProcessingUnits_min);
        if (hwinfo->nProcessingUnits_max > hwinfo->nProcessingUnits_min)
        {
            s += gmx::formatString(" - %2d", hwinfo->nProcessingUnits_max);
        }
        s += gmx::formatString("\n");
        s += gmx::formatString("  OS CPU Limit / recommended threads to start per node:   %2d",
                               hwinfo->maxThreads_min);
        if (hwinfo->maxThreads_max > hwinfo->maxThreads_min)
        {
            s += gmx::formatString(" - %2d", hwinfo->maxThreads_max);
        }
        s += gmx::formatString("\n");
        if (bGPUBinary)
        {
            s += gmx::formatString("  Compatible GPUs per node: %2d", hwinfo->ngpu_compatible_min);
            if (hwinfo->ngpu_compatible_max > hwinfo->ngpu_compatible_min)
            {
                s += gmx::formatString(" - %2d", hwinfo->ngpu_compatible_max);
            }
            s += gmx::formatString("\n");
            if (hwinfo->ngpu_compatible_tot > 0)
            {
                if (hwinfo->bIdenticalGPUs)
                {
                    s += gmx::formatString("  All nodes have identical type(s) of GPUs\n");
                }
                else
                {
                    /* This message will also appear with identical GPU types
                     * when at least one node has no GPU.
                     */
                    s += gmx::formatString(
                            "  Different nodes have different type(s) and/or order of GPUs\n");
                }
            }
        }
    }

    char host[STRLEN];
    gmx_gethostname(host, STRLEN);

#if GMX_LIB_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    s += gmx::formatString("Hardware detected on host %s (the node of MPI rank %d):\n", host, rank);
#else
    s += gmx::formatString("Hardware detected on host %s:\n", host);
#endif
    s += gmx::formatString("  CPU info:\n");

    s += gmx::formatString("    Vendor: %s\n", cpuInfo.vendorString().c_str());

    s += gmx::formatString("    Brand:  %s\n", cpuInfo.brandString().c_str());

    if (bFullCpuInfo)
    {
        s += gmx::formatString("    Family: %d   Model: %d   Stepping: %d\n",
                               cpuInfo.family(),
                               cpuInfo.model(),
                               cpuInfo.stepping());

        s += gmx::formatString("    Features:");
        for (const auto& f : cpuInfo.featureSet())
        {
            s += gmx::formatString(" %s", gmx::CpuInfo::featureString(f).c_str());
        }
        s += gmx::formatString("\n");
    }

    if (cpuInfo.feature(gmx::CpuInfo::Feature::X86_Avx512F) && cpuInfo.vendor() == gmx::CpuInfo::Vendor::Intel)
    {
        int avx512fmaunits = gmx::identifyAvx512FmaUnits();
        s += gmx::formatString("    Number of AVX-512 FMA units:");
        if (avx512fmaunits > 0)
        {
            s += gmx::formatString(" %d", avx512fmaunits);
            if (avx512fmaunits == 1)
            {
                s += gmx::formatString(" (For Intel, AVX2 is faster w/o 2 AVX-512 FMA units)");
            }
        }
        else
        {
            s += gmx::formatString(" Cannot run AVX-512 detection - assuming 2");
        }
        s += gmx::formatString("\n");
    }

    s += gmx::formatString("  Hardware topology: ");
    switch (hwTop.supportLevel())
    {
        case gmx::HardwareTopology::SupportLevel::None: s += gmx::formatString("None\n"); break;
        case gmx::HardwareTopology::SupportLevel::LogicalProcessorCount:
            s += gmx::formatString("Only logical processor count\n");
            break;
        case gmx::HardwareTopology::SupportLevel::Basic: s += gmx::formatString("Basic\n"); break;
        case gmx::HardwareTopology::SupportLevel::Full: s += gmx::formatString("Full\n"); break;
        case gmx::HardwareTopology::SupportLevel::FullWithDevices:
            s += gmx::formatString("Full, with devices\n");
            break;
    }

    if (!hwTop.isThisSystem())
    {
        s += gmx::formatString("  NOTE: Hardware topology cached or synthetic, not detected.\n");
        if (char* p = std::getenv("HWLOC_XMLFILE"))
        {
            s += gmx::formatString("        HWLOC_XMLFILE=%s\n", p);
        }
    }

    if (bFullCpuInfo)
    {
        if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Basic)
        {
            s += gmx::formatString("    Packages, cores, and logical processors:\n");
            s += gmx::formatString("    [indices refer to OS logical processors]\n");

            for (const auto& package : hwTop.machine().packages)
            {
                s += gmx::formatString("      Package %2d:", package.id);
                for (const auto& c : package.cores)
                {
                    s += gmx::formatString(" [");
                    for (const auto& pu : c.processingUnits)
                    {
                        s += gmx::formatString(" %3d", pu.osId);
                    }
                    s += gmx::formatString("]");
                }
                s += gmx::formatString("\n");
            }
        }
        s += gmx::formatString(
                "    CPU limit set by OS: %g   Recommended max number of threads: %d\n",
                hwTop.cpuLimit(),
                hwTop.maxThreads());
        if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Full)
        {
            s += gmx::formatString("    Numa nodes:\n");
            for (const auto& n : hwTop.machine().numa.nodes)
            {
                s += gmx::formatString("      Node %2d (%zu bytes mem):", n.id, n.memory);
                for (const auto& l : n.processingUnits)
                {
                    s += gmx::formatString(" %3d", hwTop.machine().logicalProcessors[l].osId);
                }
                s += gmx::formatString("\n");
            }
            s += gmx::formatString("      Latency:\n          ");
            for (std::size_t j = 0; j < hwTop.machine().numa.nodes.size(); j++)
            {
                s += gmx::formatString(" %5zu", j);
            }
            s += gmx::formatString("\n");
            for (std::size_t i = 0; i < hwTop.machine().numa.nodes.size(); i++)
            {
                s += gmx::formatString("     %5zu", i);
                for (std::size_t j = 0; j < hwTop.machine().numa.nodes.size(); j++)
                {
                    s += gmx::formatString(" %5.2f", hwTop.machine().numa.relativeLatency[i][j]);
                }
                s += gmx::formatString("\n");
            }


            s += gmx::formatString("    Caches:\n");
            for (const auto& c : hwTop.machine().caches)
            {
                s += gmx::formatString(
                        "      L%d: %zu bytes, linesize %d bytes, assoc. %d, shared %d ways\n",
                        c.level,
                        c.size,
                        c.linesize,
                        c.associativity,
                        c.shared);
            }
        }
        if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::FullWithDevices)
        {
            s += gmx::formatString("    PCI devices:\n");
            for (const auto& d : hwTop.machine().devices)
            {
                s += gmx::formatString(
                        "      %04x:%02x:%02x.%1x  Id: %04x:%04x  Class: 0x%04x  Numa: %d\n",
                        d.domain,
                        d.bus,
                        d.dev,
                        d.func,
                        d.vendorId,
                        d.deviceId,
                        d.classId,
                        d.numaNodeId);
            }
        }
    }

    if (bGPUBinary && !hwinfo->deviceInfoList.empty())
    {
        s += gmx::formatString("  GPU info:\n");
        s += gmx::formatString("    Number of GPUs detected: %d\n",
                               static_cast<int>(hwinfo->deviceInfoList.size()));
        s += sprint_gpus(hwinfo->deviceInfoList) + "\n";
    }
    return s;
}

void gmx_print_detected_hardware(FILE*                fplog,
                                 const bool           warnToStdErr,
                                 const gmx::MDLogger& mdlog,
                                 const gmx_hw_info_t* hwinfo)
{
    const gmx::CpuInfo& cpuInfo = *hwinfo->cpuInfo;

    if (fplog != nullptr)
    {
        std::string detected;

        detected = detected_hardware_string(hwinfo, TRUE);

        fprintf(fplog, "%s\n", detected.c_str());

        // Logically this should come at the end of
        // detected_hardware_string(), and can become so once that
        // function uses the MDLogger properly.
        if constexpr (bGPUBinary)
        {
            if (!hwinfo->deviceInfoList.empty())
            {
                for (const auto& deviceInfo : hwinfo->deviceInfoList)
                {
                    warnWhenDeviceNotTargeted(mdlog, *deviceInfo);
                }
            }
        }
    }

    // Do not spam stderr with all our internal information unless
    // there was something that actually went wrong; general information
    // belongs in the logfile.

    /* Check the compiled SIMD instruction set against that of the node
     * with the lowest SIMD level support (skip if SIMD detection did not work)
     */
    if (cpuInfo.supportLevel() >= gmx::CpuInfo::SupportLevel::Features)
    {
        gmx::simdCheck(cpuInfo, static_cast<gmx::SimdType>(hwinfo->simd_suggest_min), fplog, warnToStdErr);
    }

    /* For RDTSCP we only check on our local node and skip the MPI reduction, only on x86 */
    if (gmx::c_architecture == gmx::Architecture::X86)
    {
        check_use_of_rdtscp_on_this_cpu(mdlog, cpuInfo);
    }
}
