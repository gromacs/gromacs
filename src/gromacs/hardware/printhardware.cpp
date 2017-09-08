/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "printhardware.h"

#include "config.h"

#include <cstdlib>

#include <string>
#include <vector>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/simd/support.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

//! Constant used to help minimize preprocessed code
static const bool bGPUBinary     = GMX_GPU != GMX_GPU_NONE;

/*! \internal \brief
 * Returns the GPU information text, one GPU per line.
 */
static std::string sprint_gpus(const gmx_gpu_info_t &gpu_info)
{
    char                     stmp[STRLEN];
    std::vector<std::string> gpuStrings;
    for (int i = 0; i < gpu_info.n_dev; i++)
    {
        get_gpu_device_info_string(stmp, gpu_info, i);
        gpuStrings.push_back(gmx::formatString("    %s", stmp));
    }
    return gmx::joinStrings(gpuStrings, "\n");
}

/* Give a suitable fatal error or warning if the build configuration
   and runtime CPU do not match. */
static void
check_use_of_rdtscp_on_this_cpu(const gmx::MDLogger   &mdlog,
                                const gmx::CpuInfo    &cpuInfo)
{
    bool        binaryUsesRdtscp = HAVE_RDTSCP;

    const char *programName = gmx::getProgramContext().displayName();

    if (cpuInfo.supportLevel() < gmx::CpuInfo::SupportLevel::Features)
    {
        if (binaryUsesRdtscp)
        {
            GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                    "The %s executable was compiled to use the rdtscp CPU instruction. "
                    "We cannot detect the features of your current CPU, but will proceed anyway. "
                    "If you get a crash, rebuild GROMACS with the GMX_USE_RDTSCP=OFF CMake option.",
                    programName);
        }
    }
    else
    {
        bool cpuHasRdtscp = cpuInfo.feature(gmx::CpuInfo::Feature::X86_Rdtscp);

        if (!cpuHasRdtscp && binaryUsesRdtscp)
        {
            gmx_fatal(FARGS, "The %s executable was compiled to use the rdtscp CPU instruction. "
                      "However, this is not supported by the current hardware and continuing would lead to a crash. "
                      "Please rebuild GROMACS with the GMX_USE_RDTSCP=OFF CMake option.",
                      programName);
        }

        if (cpuHasRdtscp && !binaryUsesRdtscp)
        {
            GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                    "The current CPU can measure timings more accurately than the code in\n"
                    "%s was configured to use. This might affect your simulation\n"
                    "speed as accurate timings are needed for load-balancing.\n"
                    "Please consider rebuilding %s with the GMX_USE_RDTSCP=ON CMake option.",
                    programName, programName);
        }
    }
}

static std::string detected_hardware_string(const gmx_hw_info_t *hwinfo,
                                            bool                 bFullCpuInfo)
{
    std::string                  s;

    const gmx::CpuInfo          &cpuInfo = *hwinfo->cpuInfo;
    const gmx::HardwareTopology &hwTop   = *hwinfo->hardwareTopology;

    s  = gmx::formatString("\n");
    s += gmx::formatString("Running on %d node%s with total",
                           hwinfo->nphysicalnode,
                           hwinfo->nphysicalnode == 1 ? "" : "s");
    if (hwinfo->ncore_tot > 0)
    {
        s += gmx::formatString(" %d cores,", hwinfo->ncore_tot);
    }
    s += gmx::formatString(" %d logical cores", hwinfo->nhwthread_tot);
    if (hwinfo->gpu_info.bDetectGPUs)
    {
        s += gmx::formatString(", %d compatible GPU%s",
                               hwinfo->ngpu_compatible_tot,
                               hwinfo->ngpu_compatible_tot == 1 ? "" : "s");
    }
    else if (bGPUBinary)
    {
        s += gmx::formatString(" (GPU detection deactivated)");
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
        s += gmx::formatString("  Logical cores per node:   %2d", hwinfo->nhwthread_min);
        if (hwinfo->nhwthread_max > hwinfo->nhwthread_min)
        {
            s += gmx::formatString(" - %2d", hwinfo->nhwthread_max);
        }
        s += gmx::formatString("\n");
        if (bGPUBinary)
        {
            s += gmx::formatString("  Compatible GPUs per node: %2d",
                                   hwinfo->ngpu_compatible_min);
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
                    s += gmx::formatString("  Different nodes have different type(s) and/or order of GPUs\n");
                }
            }
        }
    }

#if GMX_LIB_MPI
    int  rank;
    char host[STRLEN];

    gmx_gethostname(host, STRLEN);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // TODO Use a wrapper around MPI_Get_processor_name instead.
    s += gmx::formatString("Hardware detected on host %s (the node of MPI rank %d):\n",
                           host, rank);
#else
    s += gmx::formatString("Hardware detected:\n");
#endif
    s += gmx::formatString("  CPU info:\n");

    s += gmx::formatString("    Vendor: %s\n", cpuInfo.vendorString().c_str());

    s += gmx::formatString("    Brand:  %s\n", cpuInfo.brandString().c_str());

    if (bFullCpuInfo)
    {
        s += gmx::formatString("    Family: %d   Model: %d   Stepping: %d\n",
                               cpuInfo.family(), cpuInfo.model(), cpuInfo.stepping());

        s += gmx::formatString("    Features:");
        for (auto &f : cpuInfo.featureSet())
        {
            s += gmx::formatString(" %s", cpuInfo.featureString(f).c_str());;
        }
        s += gmx::formatString("\n");
    }

    s += gmx::formatString("    SIMD instructions most likely to fit this hardware: %s",
                           gmx::simdString(static_cast<gmx::SimdType>(hwinfo->simd_suggest_min)).c_str());

    if (hwinfo->simd_suggest_max > hwinfo->simd_suggest_min)
    {
        s += gmx::formatString(" - %s", gmx::simdString(static_cast<gmx::SimdType>(hwinfo->simd_suggest_max)).c_str());
    }
    s += gmx::formatString("\n");

    s += gmx::formatString("    SIMD instructions selected at GROMACS compile time: %s\n",
                           gmx::simdString(gmx::simdCompiled()).c_str());

    s += gmx::formatString("\n");

    s += gmx::formatString("  Hardware topology: ");
    switch (hwTop.supportLevel())
    {
        case gmx::HardwareTopology::SupportLevel::None:
            s += gmx::formatString("None\n");
            break;
        case gmx::HardwareTopology::SupportLevel::LogicalProcessorCount:
            s += gmx::formatString("Only logical processor count\n");
            break;
        case gmx::HardwareTopology::SupportLevel::Basic:
            s += gmx::formatString("Basic\n");
            break;
        case gmx::HardwareTopology::SupportLevel::Full:
            s += gmx::formatString("Full\n");
            break;
        case gmx::HardwareTopology::SupportLevel::FullWithDevices:
            s += gmx::formatString("Full, with devices\n");
            break;
    }

    if (!hwTop.isThisSystem())
    {
        s += gmx::formatString("  NOTE: Hardware topology cached or synthetic, not detected.\n");
        if (char *p = std::getenv("HWLOC_XMLFILE"))
        {
            s += gmx::formatString("        HWLOC_XMLFILE=%s\n", p);
        }
    }

    if (bFullCpuInfo)
    {
        if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Basic)
        {
            s += gmx::formatString("    Sockets, cores, and logical processors:\n");

            for (auto &socket : hwTop.machine().sockets)
            {
                s += gmx::formatString("      Socket %2d:", socket.id);
                for (auto &c : socket.cores)
                {
                    s += gmx::formatString(" [");
                    for (auto &t : c.hwThreads)
                    {
                        s += gmx::formatString(" %3d", t.logicalProcessorId);
                    }
                    s += gmx::formatString("]");
                }
                s += gmx::formatString("\n");
            }
        }
        if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Full)
        {
            s += gmx::formatString("    Numa nodes:\n");
            for (auto &n : hwTop.machine().numa.nodes)
            {
                s += gmx::formatString("      Node %2d (%" GMX_PRIu64 " bytes mem):", n.id, n.memory);
                for (auto &l : n.logicalProcessorId)
                {
                    s += gmx::formatString(" %3d", l);
                }
                s += gmx::formatString("\n");
            }
            s += gmx::formatString("      Latency:\n          ");
            for (std::size_t j = 0; j < hwTop.machine().numa.nodes.size(); j++)
            {
                s += gmx::formatString(" %5d", j);
            }
            s += gmx::formatString("\n");
            for (std::size_t i = 0; i < hwTop.machine().numa.nodes.size(); i++)
            {
                s += gmx::formatString("     %5d", i);
                for (std::size_t j = 0; j < hwTop.machine().numa.nodes.size(); j++)
                {
                    s += gmx::formatString(" %5.2f", hwTop.machine().numa.relativeLatency[i][j]);
                }
                s += gmx::formatString("\n");
            }


            s += gmx::formatString("    Caches:\n");
            for (auto &c : hwTop.machine().caches)
            {
                s += gmx::formatString("      L%d: %" GMX_PRIu64 " bytes, linesize %d bytes, assoc. %d, shared %d ways\n",
                                       c.level, c.size, c.linesize, c.associativity, c.shared);
            }
        }
        if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::FullWithDevices)
        {
            s += gmx::formatString("    PCI devices:\n");
            for (auto &d : hwTop.machine().devices)
            {
                s += gmx::formatString("      %04x:%02x:%02x.%1x  Id: %04x:%04x  Class: 0x%04x  Numa: %d\n",
                                       d.domain, d.bus, d.dev, d.func, d.vendorId, d.deviceId, d.classId, d.numaNodeId);
            }
        }
    }

    if (bGPUBinary && (hwinfo->ngpu_compatible_tot > 0 ||
                       hwinfo->gpu_info.n_dev > 0))
    {
        s += gmx::formatString("  GPU info:\n");
        s += gmx::formatString("    Number of GPUs detected: %d\n",
                               hwinfo->gpu_info.n_dev);
        if (hwinfo->gpu_info.n_dev > 0)
        {
            s += sprint_gpus(hwinfo->gpu_info) + "\n";
        }
    }
    return s;
}

void gmx_print_detected_hardware(FILE *fplog, const t_commrec *cr,
                                 const gmx::MDLogger &mdlog,
                                 const gmx_hw_info_t *hwinfo)
{
    const gmx::CpuInfo &cpuInfo = *hwinfo->cpuInfo;

    if (fplog != nullptr)
    {
        std::string detected;

        detected = detected_hardware_string(hwinfo, TRUE);

        fprintf(fplog, "%s\n", detected.c_str());
    }

    if (MULTIMASTER(cr))
    {
        std::string detected;

        detected = detected_hardware_string(hwinfo, FALSE);

        fprintf(stderr, "%s\n", detected.c_str());
    }

    /* Check the compiled SIMD instruction set against that of the node
     * with the lowest SIMD level support (skip if SIMD detection did not work)
     */
    if (cpuInfo.supportLevel() >= gmx::CpuInfo::SupportLevel::Features)
    {
        gmx::simdCheck(static_cast<gmx::SimdType>(hwinfo->simd_suggest_min), fplog, MULTIMASTER(cr));
    }

    /* For RDTSCP we only check on our local node and skip the MPI reduction */
    check_use_of_rdtscp_on_this_cpu(mdlog, cpuInfo);
}
