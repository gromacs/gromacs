/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

#include <cerrno>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <chrono>
#include <string>
#include <thread>
#include <vector>

#include "thread_mpi/threads.h"

#include "gromacs/gmxlib/md_logging.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/simd/support.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#ifdef HAVE_UNISTD_H
#    include <unistd.h>       // sysconf()
#endif

//! Convenience macro to help us avoid ifdefs each time we use sysconf
#if !defined(_SC_NPROCESSORS_ONLN) && defined(_SC_NPROC_ONLN)
#    define _SC_NPROCESSORS_ONLN _SC_NPROC_ONLN
#endif

//! Convenience macro to help us avoid ifdefs each time we use sysconf
#if !defined(_SC_NPROCESSORS_CONF) && defined(_SC_NPROC_CONF)
#    define _SC_NPROCESSORS_CONF _SC_NPROC_CONF
#endif

#if defined (__i386__) || defined (__x86_64__) || defined (_M_IX86) || defined (_M_X64)
//! Constant used to help minimize preprocessed code
static const bool isX86 = true;
#else
//! Constant used to help minimize preprocessed code
static const bool isX86 = false;
#endif

#if defined __powerpc__ || defined __ppc__ || defined __PPC__
static const bool isPowerPC = true;
#else
static const bool isPowerPC = false;
#endif

//! Constant used to help minimize preprocessed code
static const bool bGPUBinary     = GMX_GPU != GMX_GPU_NONE;

/* Note that some of the following arrays must match the "GPU support
 * enumeration" in src/config.h.cmakein, so that GMX_GPU looks up an
 * array entry. */

/* Both CUDA and OpenCL (on the supported/tested platforms) supports
 * GPU device sharing.
 */
static const bool gpuSharingSupport[] = { false, true, true };
static const bool bGpuSharingSupported = gpuSharingSupport[GMX_GPU];

/* Both CUDA and OpenCL (on the tested/supported platforms) supports everything.
 */
static const bool multiGpuSupport[] = {
    false, true, true
};
static const bool bMultiGpuPerNodeSupported = multiGpuSupport[GMX_GPU];

/* Names of the GPU detection/check results (see e_gpu_detect_res_t in hw_info.h). */
const char * const gpu_detect_res_str[egpuNR] =
{
    "compatible", "inexistent", "incompatible", "insane"
};

static const char * invalid_gpuid_hint =
    "A delimiter-free sequence of valid numeric IDs of available GPUs is expected.";

/* The globally shared hwinfo structure. */
static gmx_hw_info_t      *hwinfo_g;
/* A reference counter for the hwinfo structure */
static int                 n_hwinfo = 0;
/* A lock to protect the hwinfo structure */
static tMPI_Thread_mutex_t hw_info_lock = TMPI_THREAD_MUTEX_INITIALIZER;

#define HOSTNAMELEN 80

/* FW decl. */
static void set_gpu_ids(gmx_gpu_opt_t *gpu_opt, int nrank, int rank);
static int gmx_count_gpu_dev_unique(const gmx_gpu_info_t *gpu_info,
                                    const gmx_gpu_opt_t  *gpu_opt);

gmx_bool gmx_multiple_gpu_per_node_supported()
{
    return bMultiGpuPerNodeSupported;
}

gmx_bool gmx_gpu_sharing_supported()
{
    return bGpuSharingSupported;
}

static void sprint_gpus(char *sbuf, const gmx_gpu_info_t *gpu_info)
{
    int      i, ndev;
    char     stmp[STRLEN];

    ndev = gpu_info->n_dev;

    sbuf[0] = '\0';
    for (i = 0; i < ndev; i++)
    {
        get_gpu_device_info_string(stmp, gpu_info, i);
        strcat(sbuf, "    ");
        strcat(sbuf, stmp);
        if (i < ndev - 1)
        {
            strcat(sbuf, "\n");
        }
    }
}

static void print_gpu_detection_stats(FILE                 *fplog,
                                      const gmx_gpu_info_t *gpu_info,
                                      const t_commrec      *cr)
{
    char onhost[HOSTNAMELEN+10], stmp[STRLEN];
    int  ngpu;

    if (!gpu_info->bDetectGPUs)
    {
        /* We skipped the detection, so don't print detection stats */
        return;
    }

    ngpu = gpu_info->n_dev;

#if GMX_LIB_MPI
    /* We only print the detection on one, of possibly multiple, nodes */
    std::strncpy(onhost, " on host ", 10);
    gmx_gethostname(onhost + 9, HOSTNAMELEN);
#else
    /* We detect all relevant GPUs */
    std::strncpy(onhost, "", 1);
#endif

    if (ngpu > 0)
    {
        sprint_gpus(stmp, gpu_info);
        md_print_warn(cr, fplog, "%d GPU%s detected%s:\n%s\n",
                      ngpu, (ngpu > 1) ? "s" : "", onhost, stmp);
    }
    else
    {
        md_print_warn(cr, fplog, "No GPUs detected%s\n", onhost);
    }
}

/*! \brief Helper function for reporting GPU usage information
 * in the mdrun log file
 *
 * \param[in] gpu_info       Pointer to per-node GPU info struct
 * \param[in] gpu_opt        Pointer to per-node GPU options struct
 * \param[in] numPpRanks     Number of PP ranks per node
 * \param[in] bPrintHostName Print the hostname in the usage information
 * \return                   String to write to the log file
 * \throws                   std::bad_alloc if out of memory */
static std::string
makeGpuUsageReport(const gmx_gpu_info_t *gpu_info,
                   const gmx_gpu_opt_t  *gpu_opt,
                   size_t                numPpRanks,
                   bool                  bPrintHostName)
{
    int  ngpu_use  = gpu_opt->n_dev_use;
    int  ngpu_comp = gpu_info->n_dev_compatible;
    char host[HOSTNAMELEN];

    if (bPrintHostName)
    {
        gmx_gethostname(host, HOSTNAMELEN);
    }

    /* Issue a note if GPUs are available but not used */
    if (ngpu_comp > 0 && ngpu_use < 1)
    {
        return gmx::formatString("%d compatible GPU%s detected in the system, but none will be used.\n"
                                 "Consider trying GPU acceleration with the Verlet scheme!\n",
                                 ngpu_comp, (ngpu_comp > 1) ? "s" : "");
    }

    std::string output;
    if (!gpu_opt->bUserSet)
    {
        // gpu_opt->dev_compatible is only populated during auto-selection
        std::string gpuIdsString =
            formatAndJoin(gmx::constArrayRefFromArray(gpu_opt->dev_compatible,
                                                      gpu_opt->n_dev_compatible),
                          ",", gmx::StringFormatter("%d"));
        bool bPluralGpus = gpu_opt->n_dev_compatible > 1;

        if (bPrintHostName)
        {
            output += gmx::formatString("On host %s ", host);
        }
        output += gmx::formatString("%d compatible GPU%s %s present, with ID%s %s\n",
                                    gpu_opt->n_dev_compatible,
                                    bPluralGpus ? "s" : "",
                                    bPluralGpus ? "are" : "is",
                                    bPluralGpus ? "s" : "",
                                    gpuIdsString.c_str());
    }

    {
        std::vector<int> gpuIdsInUse;
        for (int i = 0; i < ngpu_use; i++)
        {
            gpuIdsInUse.push_back(get_gpu_device_id(gpu_info, gpu_opt, i));
        }
        std::string gpuIdsString =
            formatAndJoin(gpuIdsInUse, ",", gmx::StringFormatter("%d"));
        int         numGpusInUse = gmx_count_gpu_dev_unique(gpu_info, gpu_opt);
        bool        bPluralGpus  = numGpusInUse > 1;

        if (bPrintHostName)
        {
            output += gmx::formatString("On host %s ", host);
        }
        output += gmx::formatString("%d GPU%s %sselected for this run.\n"
                                    "Mapping of GPU ID%s to the %d PP rank%s in this node: %s\n",
                                    numGpusInUse, bPluralGpus ? "s" : "",
                                    gpu_opt->bUserSet ? "user-" : "auto-",
                                    bPluralGpus ? "s" : "",
                                    numPpRanks,
                                    (numPpRanks > 1) ? "s" : "",
                                    gpuIdsString.c_str());
    }

    return output;
}

/* Give a suitable fatal error or warning if the build configuration
   and runtime CPU do not match. */
static void
check_use_of_rdtscp_on_this_cpu(FILE                  *fplog,
                                const t_commrec       *cr,
                                const gmx::CpuInfo    &cpuInfo)
{
#ifdef HAVE_RDTSCP
    bool binaryUsesRdtscp = TRUE;
#else
    bool binaryUsesRdtscp = FALSE;
#endif

    const char *programName = gmx::getProgramContext().displayName();

    if (cpuInfo.supportLevel() < gmx::CpuInfo::SupportLevel::Features)
    {
        if (binaryUsesRdtscp)
        {
            md_print_warn(cr, fplog, "The %s executable was compiled to use the rdtscp CPU instruction. "
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
            md_print_warn(cr, fplog, "The current CPU can measure timings more accurately than the code in\n"
                          "%s was configured to use. This might affect your simulation\n"
                          "speed as accurate timings are needed for load-balancing.\n"
                          "Please consider rebuilding %s with the GMX_USE_RDTSCP=ON CMake option.\n",
                          programName, programName);
        }
    }
}

void gmx_check_hw_runconf_consistency(FILE                *fplog,
                                      const gmx_hw_info_t *hwinfo,
                                      const t_commrec     *cr,
                                      const gmx_hw_opt_t  *hw_opt,
                                      gmx_bool             bUseGPU)
{
    int      npppn;
    char     th_or_proc[STRLEN], th_or_proc_plural[STRLEN], pernode[STRLEN];
    gmx_bool btMPI, bMPI, bNthreadsAuto, bEmulateGPU;

    GMX_RELEASE_ASSERT(hwinfo, "hwinfo must be a non-NULL pointer");
    GMX_RELEASE_ASSERT(cr, "cr must be a non-NULL pointer");

    /* Below we only do consistency checks for PP and GPUs,
     * this is irrelevant for PME only nodes, so in that case we return
     * here.
     */
    if (!(cr->duty & DUTY_PP))
    {
        return;
    }

#if GMX_THREAD_MPI
    bMPI          = FALSE;
    btMPI         = TRUE;
    bNthreadsAuto = (hw_opt->nthreads_tmpi < 1);
#elif GMX_LIB_MPI
    bMPI          = TRUE;
    btMPI         = FALSE;
    bNthreadsAuto = FALSE;
#else
    bMPI          = FALSE;
    btMPI         = FALSE;
    bNthreadsAuto = FALSE;
#endif

    /* GPU emulation detection is done later, but we need here as well
     * -- uncool, but there's no elegant workaround */
    bEmulateGPU       = (getenv("GMX_EMULATE_GPU") != NULL);

    if (hwinfo->gpu_info.n_dev_compatible > 0)
    {
        std::string gpuUseageReport;
        try
        {
            gpuUseageReport = makeGpuUsageReport(&hwinfo->gpu_info,
                                                 &hw_opt->gpu_opt,
                                                 cr->nrank_pp_intranode,
                                                 bMPI && cr->nnodes > 1);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

        /* NOTE: this print is only for and on one physical node */
        md_print_info(cr, fplog, "%s\n", gpuUseageReport.c_str());
    }

    /* Need to ensure that we have enough GPUs:
     * - need one GPU per PP node
     * - no GPU oversubscription with tMPI
     * */
    /* number of PP processes per node */
    npppn = cr->nrank_pp_intranode;

    pernode[0]           = '\0';
    th_or_proc_plural[0] = '\0';
    if (btMPI)
    {
        sprintf(th_or_proc, "thread-MPI thread");
        if (npppn > 1)
        {
            sprintf(th_or_proc_plural, "s");
        }
    }
    else if (bMPI)
    {
        sprintf(th_or_proc, "MPI process");
        if (npppn > 1)
        {
            sprintf(th_or_proc_plural, "es");
        }
        sprintf(pernode, " per node");
    }
    else
    {
        /* neither MPI nor tMPI */
        sprintf(th_or_proc, "process");
    }

    if (bUseGPU && hwinfo->gpu_info.n_dev_compatible > 0 &&
        !bEmulateGPU)
    {
        int  ngpu_comp, ngpu_use;
        char gpu_comp_plural[2], gpu_use_plural[2];

        ngpu_comp = hwinfo->gpu_info.n_dev_compatible;
        ngpu_use  = hw_opt->gpu_opt.n_dev_use;

        sprintf(gpu_comp_plural, "%s", (ngpu_comp > 1) ? "s" : "");
        sprintf(gpu_use_plural,  "%s", (ngpu_use > 1) ? "s" : "");

        const char *programName = gmx::getProgramContext().displayName();

        /* number of tMPI threads auto-adjusted */
        if (btMPI && bNthreadsAuto)
        {
            if (hw_opt->gpu_opt.bUserSet && npppn < ngpu_use)
            {
                /* The user manually provided more GPUs than threads we
                   could automatically start. */
                gmx_fatal(FARGS,
                          "%d GPU%s provided, but only %d PP thread-MPI thread%s coud be started.\n"
                          "%s requires one PP tread-MPI thread per GPU; use fewer GPUs.",
                          ngpu_use, gpu_use_plural,
                          npppn, th_or_proc_plural,
                          programName);
            }

            if (!hw_opt->gpu_opt.bUserSet && npppn < ngpu_comp)
            {
                /* There are more GPUs than tMPI threads; we have
                   limited the number GPUs used. */
                md_print_warn(cr, fplog,
                              "NOTE: %d GPU%s were detected, but only %d PP thread-MPI thread%s can be started.\n"
                              "      %s can use one GPU per PP tread-MPI thread, so only %d GPU%s will be used.\n",
                              ngpu_comp, gpu_comp_plural,
                              npppn, th_or_proc_plural,
                              programName, npppn,
                              npppn > 1 ? "s" : "");
            }
        }

        if (hw_opt->gpu_opt.bUserSet)
        {
            if (ngpu_use != npppn)
            {
                gmx_fatal(FARGS,
                          "Incorrect launch configuration: mismatching number of PP %s%s and GPUs%s.\n"
                          "%s was started with %d PP %s%s%s, but you provided %d GPU%s.",
                          th_or_proc, btMPI ? "s" : "es", pernode,
                          programName, npppn, th_or_proc,
                          th_or_proc_plural, pernode,
                          ngpu_use, gpu_use_plural);
            }
        }
        else
        {
            /* TODO Should we have a gpu_opt->n_dev_supported field? */
            if (ngpu_comp > npppn && gmx_multiple_gpu_per_node_supported())
            {
                md_print_warn(cr, fplog,
                              "NOTE: potentially sub-optimal launch configuration, %s started with less\n"
                              "      PP %s%s%s than GPU%s available.\n"
                              "      Each PP %s can use only one GPU, %d GPU%s%s will be used.\n",
                              programName, th_or_proc,
                              th_or_proc_plural, pernode, gpu_comp_plural,
                              th_or_proc, npppn, gpu_use_plural, pernode);
            }

            if (ngpu_use != npppn)
            {
                /* Avoid duplicate error messages.
                 * Unfortunately we can only do this at the physical node
                 * level, since the hardware setup and MPI process count
                 * might differ between physical nodes.
                 */
                if (cr->rank_pp_intranode == 0)
                {
                    std::string reasonForLimit;
                    if (ngpu_comp > 1 &&
                        ngpu_use == 1 &&
                        !gmx_multiple_gpu_per_node_supported())
                    {
                        reasonForLimit  = "can be used by ";
                        reasonForLimit += getGpuImplementationString();
                        reasonForLimit += " in GROMACS";
                    }
                    else
                    {
                        reasonForLimit = "was detected";
                    }
                    gmx_fatal(FARGS,
                              "Incorrect launch configuration: mismatching number of PP %s%s and GPUs%s.\n"
                              "%s was started with %d PP %s%s%s, but only %d GPU%s %s.",
                              th_or_proc, btMPI ? "s" : "es", pernode,
                              programName, npppn, th_or_proc,
                              th_or_proc_plural, pernode,
                              ngpu_use, gpu_use_plural, reasonForLimit.c_str());
                }
            }
        }

        {
            int      same_count;

            same_count = gmx_count_gpu_dev_shared(&hw_opt->gpu_opt);

            if (same_count > 0)
            {
                md_print_info(cr, fplog,
                              "NOTE: You assigned %s to multiple %s%s.\n",
                              same_count > 1 ? "GPUs" : "a GPU", th_or_proc, btMPI ? "s" : "es");
            }
        }
    }

#if GMX_MPI
    if (PAR(cr))
    {
        /* Avoid other ranks to continue after
           inconsistency */
        MPI_Barrier(cr->mpi_comm_mygroup);
    }
#endif

}

/* Return 0 if none of the GPU (per node) are shared among PP ranks.
 *
 * Sharing GPUs among multiple PP ranks is possible when the user passes
 * GPU IDs. Here we check for sharing and return a non-zero value when
 * this is detected. Note that the return value represents the number of
 * PP rank pairs that share a device.
 */
int gmx_count_gpu_dev_shared(const gmx_gpu_opt_t *gpu_opt)
{
    int      same_count    = 0;
    int      ngpu          = gpu_opt->n_dev_use;

    if (gpu_opt->bUserSet)
    {
        int      i, j;

        for (i = 0; i < ngpu - 1; i++)
        {
            for (j = i + 1; j < ngpu; j++)
            {
                same_count      += (gpu_opt->dev_use[i] ==
                                    gpu_opt->dev_use[j]);
            }
        }
    }

    return same_count;
}

/* Count and return the number of unique GPUs (per node) selected.
 *
 * As sharing GPUs among multiple PP ranks is possible when the user passes
 * GPU IDs, the number of GPUs user (per node) can be different from the
 * number of GPU IDs selected.
 */
static int gmx_count_gpu_dev_unique(const gmx_gpu_info_t *gpu_info,
                                    const gmx_gpu_opt_t  *gpu_opt)
{
    int  i, uniq_count, ngpu;
    int *uniq_ids;

    GMX_RELEASE_ASSERT(gpu_info, "gpu_info must be a non-NULL pointer");
    GMX_RELEASE_ASSERT(gpu_opt, "gpu_opt must be a non-NULL pointer");

    ngpu = gpu_info->n_dev;

    uniq_count  = 0;

    snew(uniq_ids, ngpu);

    /* Each element in uniq_ids will be set to 0 or 1. The n-th element set
     * to 1 indicates that the respective GPU was selected to be used. */
    for (i = 0; i < gpu_opt->n_dev_use; i++)
    {
        int device_id;

        device_id           = gmx_gpu_sharing_supported() ? get_gpu_device_id(gpu_info, gpu_opt, i) : i;
        uniq_ids[device_id] = 1;
    }
    /* Count the devices used. */
    for (i = 0; i < ngpu; i++)
    {
        uniq_count += uniq_ids[i];
    }

    sfree(uniq_ids);

    return uniq_count;
}

static void gmx_detect_gpus(FILE *fplog, const t_commrec *cr)
{
#if GMX_LIB_MPI
    int              rank_world;
    MPI_Comm         physicalnode_comm;
#endif
    int              rank_local;

    /* Under certain circumstances MPI ranks on the same physical node
     * can not simultaneously access the same GPU(s). Therefore we run
     * the detection only on one MPI rank per node and broadcast the info.
     * Note that with thread-MPI only a single thread runs this code.
     *
     * NOTE: We can't broadcast gpu_info with OpenCL as the device and platform
     * ID stored in the structure are unique for each rank (even if a device
     * is shared by multiple ranks).
     *
     * TODO: We should also do CPU hardware detection only once on each
     * physical node and broadcast it, instead of do it on every MPI rank.
     */
#if GMX_LIB_MPI
    /* A split of MPI_COMM_WORLD over physical nodes is only required here,
     * so we create and destroy it locally.
     */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
    MPI_Comm_split(MPI_COMM_WORLD, gmx_physicalnode_id_hash(),
                   rank_world, &physicalnode_comm);
    MPI_Comm_rank(physicalnode_comm, &rank_local);
#else
    /* Here there should be only one process, check this */
    GMX_RELEASE_ASSERT(cr->nnodes == 1 && cr->sim_nodeid == 0, "Only a single (master) process should execute here");

    rank_local = 0;
#endif

    /*  With CUDA detect only on one rank per host, with OpenCL need do
     *  the detection on all PP ranks */
    bool isOpenclPpRank = ((GMX_GPU == GMX_GPU_OPENCL) && (cr->duty & DUTY_PP));

    if (rank_local == 0 || isOpenclPpRank)
    {
        char detection_error[STRLEN] = "", sbuf[STRLEN];

        if (detect_gpus(&hwinfo_g->gpu_info, detection_error) != 0)
        {
            if (detection_error[0] != '\0')
            {
                sprintf(sbuf, ":\n      %s\n", detection_error);
            }
            else
            {
                sprintf(sbuf, ".");
            }
            md_print_warn(cr, fplog,
                          "NOTE: Error occurred during GPU detection%s"
                          "      Can not use GPU acceleration, will fall back to CPU kernels.\n",
                          sbuf);
        }
    }

#if GMX_LIB_MPI
    if (!isOpenclPpRank)
    {
        /* Broadcast the GPU info to the other ranks within this node */
        MPI_Bcast(&hwinfo_g->gpu_info.n_dev, 1, MPI_INT, 0, physicalnode_comm);

        if (hwinfo_g->gpu_info.n_dev > 0)
        {
            int dev_size;

            dev_size = hwinfo_g->gpu_info.n_dev*sizeof_gpu_dev_info();

            if (rank_local > 0)
            {
                hwinfo_g->gpu_info.gpu_dev =
                    (struct gmx_device_info_t *)malloc(dev_size);
            }
            MPI_Bcast(hwinfo_g->gpu_info.gpu_dev, dev_size, MPI_BYTE,
                      0, physicalnode_comm);
            MPI_Bcast(&hwinfo_g->gpu_info.n_dev_compatible, 1, MPI_INT,
                      0, physicalnode_comm);
        }
    }

    MPI_Comm_free(&physicalnode_comm);
#endif
}

static void gmx_collect_hardware_mpi(const gmx::CpuInfo &cpuInfo)
{
    const int ncore = hwinfo_g->hardwareTopology->numberOfCores();
#if GMX_LIB_MPI
    int       rank_id;
    int       nrank, rank, nhwthread, ngpu, i;
    int       gpu_hash;
    int      *buf, *all;

    rank_id   = gmx_physicalnode_id_hash();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    nhwthread = hwinfo_g->nthreads_hw_avail;
    ngpu      = hwinfo_g->gpu_info.n_dev_compatible;
    /* Create a unique hash of the GPU type(s) in this node */
    gpu_hash  = 0;
    /* Here it might be better to only loop over the compatible GPU, but we
     * don't have that information available and it would also require
     * removing the device ID from the device info string.
     */
    for (i = 0; i < hwinfo_g->gpu_info.n_dev; i++)
    {
        char stmp[STRLEN];

        /* Since the device ID is incorporated in the hash, the order of
         * the GPUs affects the hash. Also two identical GPUs won't give
         * a gpu_hash of zero after XORing.
         */
        get_gpu_device_info_string(stmp, &hwinfo_g->gpu_info, i);
        gpu_hash ^= gmx_string_fullhash_func(stmp, gmx_string_hash_init);
    }

    snew(buf, nrank);
    snew(all, nrank);
    buf[rank] = rank_id;

    MPI_Allreduce(buf, all, nrank, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    gmx_bool bFound;
    int      nnode0, ncore0, nhwthread0, ngpu0, r;

    bFound     = FALSE;
    ncore0     = 0;
    nnode0     = 0;
    nhwthread0 = 0;
    ngpu0      = 0;
    for (r = 0; r < nrank; r++)
    {
        if (all[r] == rank_id)
        {
            if (!bFound && r == rank)
            {
                /* We are the first rank in this physical node */
                nnode0     = 1;
                ncore0     = ncore;
                nhwthread0 = nhwthread;
                ngpu0      = ngpu;
            }
            bFound = TRUE;
        }
    }

    sfree(buf);
    sfree(all);

    int sum[4], maxmin[10];

    {
        int buf[4];

        /* Sum values from only intra-rank 0 so we get the sum over all nodes */
        buf[0] = nnode0;
        buf[1] = ncore0;
        buf[2] = nhwthread0;
        buf[3] = ngpu0;

        MPI_Allreduce(buf, sum, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    {
        int buf[10];

        /* Store + and - values for all ranks,
         * so we can get max+min with one MPI call.
         */
        buf[0] = ncore;
        buf[1] = nhwthread;
        buf[2] = ngpu;
        buf[3] = static_cast<int>(gmx::simdSuggested(cpuInfo));
        buf[4] = gpu_hash;
        buf[5] = -buf[0];
        buf[6] = -buf[1];
        buf[7] = -buf[2];
        buf[8] = -buf[3];
        buf[9] = -buf[4];

        MPI_Allreduce(buf, maxmin, 10, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }

    hwinfo_g->nphysicalnode       = sum[0];
    hwinfo_g->ncore_tot           = sum[1];
    hwinfo_g->ncore_min           = -maxmin[5];
    hwinfo_g->ncore_max           = maxmin[0];
    hwinfo_g->nhwthread_tot       = sum[2];
    hwinfo_g->nhwthread_min       = -maxmin[6];
    hwinfo_g->nhwthread_max       = maxmin[1];
    hwinfo_g->ngpu_compatible_tot = sum[3];
    hwinfo_g->ngpu_compatible_min = -maxmin[7];
    hwinfo_g->ngpu_compatible_max = maxmin[2];
    hwinfo_g->simd_suggest_min    = -maxmin[8];
    hwinfo_g->simd_suggest_max    = maxmin[3];
    hwinfo_g->bIdenticalGPUs      = (maxmin[4] == -maxmin[9]);
#else
    /* All ranks use the same pointer, protected by a mutex in the caller */
    hwinfo_g->nphysicalnode       = 1;
    hwinfo_g->ncore_tot           = ncore;
    hwinfo_g->ncore_min           = ncore;
    hwinfo_g->ncore_max           = ncore;
    hwinfo_g->nhwthread_tot       = hwinfo_g->nthreads_hw_avail;
    hwinfo_g->nhwthread_min       = hwinfo_g->nthreads_hw_avail;
    hwinfo_g->nhwthread_max       = hwinfo_g->nthreads_hw_avail;
    hwinfo_g->ngpu_compatible_tot = hwinfo_g->gpu_info.n_dev_compatible;
    hwinfo_g->ngpu_compatible_min = hwinfo_g->gpu_info.n_dev_compatible;
    hwinfo_g->ngpu_compatible_max = hwinfo_g->gpu_info.n_dev_compatible;
    hwinfo_g->simd_suggest_min    = static_cast<int>(simdSuggested(cpuInfo));
    hwinfo_g->simd_suggest_max    = static_cast<int>(simdSuggested(cpuInfo));
    hwinfo_g->bIdenticalGPUs      = TRUE;
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
static void
spinUpCore() noexcept
{
#if defined(HAVE_SYSCONF) && defined(_SC_NPROCESSORS_CONF) && defined(_SC_NPROCESSORS_ONLN)
    // steady_clock is better than system_clock, but unsupported in gcc-4.6.4.
    // For release-2017 we can retire gcc-4.6 support and move to steady_clock.
    float dummy           = 0.1;
    int   countConfigured = sysconf(_SC_NPROCESSORS_CONF);    // noexcept
    auto  start           = std::chrono::system_clock::now(); // noexcept

    while (sysconf(_SC_NPROCESSORS_ONLN) < countConfigured &&
           std::chrono::system_clock::now() - start < std::chrono::seconds(2))
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
static void
hardwareTopologyPrepareDetection()
{
#if defined(HAVE_SYSCONF) && defined(_SC_NPROCESSORS_CONF) && \
    (defined(THREAD_PTHREADS) || defined(THREAD_WINDOWS))

    // Modify this conditional when/if x86 or PowerPC starts to sleep some cores
    if (!isX86 && !isPowerPC)
    {
        int                      countConfigured  = sysconf(_SC_NPROCESSORS_CONF);
        std::vector<std::thread> workThreads(countConfigured);

        for (auto &t : workThreads)
        {
            t = std::thread(spinUpCore);
        }

        for (auto &t : workThreads)
        {
            t.join();
        }
    }
#endif
}

/*! \brief Sanity check hardware topology and optionally print some notes to log
 *
 *  \param fplog            Log file pointer. This can be NULL, but the then routine
 *                          will not do anything.
 *  \param hardwareTopology Reference to hardwareTopology object.
 */
static void
hardwareTopologyDoubleCheckDetection(FILE gmx_unused                          *fplog,
                                     const gmx::HardwareTopology gmx_unused   &hardwareTopology)
{
#if defined HAVE_SYSCONF && defined(_SC_NPROCESSORS_CONF)
    if (fplog == NULL ||
        hardwareTopology.supportLevel() < gmx::HardwareTopology::SupportLevel::LogicalProcessorCount)
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
        fprintf(fplog, "Note: %d CPUs configured, but only %d were detected to be online.\n", countConfigured, countFromDetection);

        if (isX86 && countConfigured == 2*countFromDetection)
        {
            fprintf(fplog, "      X86 Hyperthreading is likely disabled; enable it for better performance.\n");
        }
        // For PowerPC (likely Power8) it is possible to set SMT to either 2,4, or 8-way hardware threads.
        // We only warn if it is completely disabled since default performance drops with SMT8.
        if (isPowerPC && countConfigured == 8*countFromDetection)
        {
            fprintf(fplog, "      PowerPC SMT is likely disabled; enable SMT2/SMT4 for better performance.\n");
        }
    }
#endif
}


gmx_hw_info_t *gmx_detect_hardware(FILE *fplog, const t_commrec *cr,
                                   gmx_bool bDetectGPUs)
{
    int ret;

    /* make sure no one else is doing the same thing */
    ret = tMPI_Thread_mutex_lock(&hw_info_lock);
    if (ret != 0)
    {
        gmx_fatal(FARGS, "Error locking hwinfo mutex: %s", strerror(errno));
    }

    /* only initialize the hwinfo structure if it is not already initalized */
    if (n_hwinfo == 0)
    {
        snew(hwinfo_g, 1);

        hwinfo_g->cpuInfo             = new gmx::CpuInfo(gmx::CpuInfo::detect());

        hardwareTopologyPrepareDetection();
        hwinfo_g->hardwareTopology    = new gmx::HardwareTopology(gmx::HardwareTopology::detect());

        // If we detected the topology on this system, double-check that it makes sense
        if (hwinfo_g->hardwareTopology->isThisSystem())
        {
            hardwareTopologyDoubleCheckDetection(fplog, *(hwinfo_g->hardwareTopology));
        }

        // TODO: Get rid of this altogether.
        hwinfo_g->nthreads_hw_avail = hwinfo_g->hardwareTopology->machine().logicalProcessorCount;

        /* detect GPUs */
        hwinfo_g->gpu_info.n_dev            = 0;
        hwinfo_g->gpu_info.n_dev_compatible = 0;
        hwinfo_g->gpu_info.gpu_dev          = NULL;

        /* Run the detection if the binary was compiled with GPU support
         * and we requested detection.
         */
        hwinfo_g->gpu_info.bDetectGPUs =
            (bGPUBinary && bDetectGPUs &&
             getenv("GMX_DISABLE_GPU_DETECTION") == NULL);
        if (hwinfo_g->gpu_info.bDetectGPUs)
        {
            gmx_detect_gpus(fplog, cr);
        }

        gmx_collect_hardware_mpi(*hwinfo_g->cpuInfo);
    }
    /* increase the reference counter */
    n_hwinfo++;

    ret = tMPI_Thread_mutex_unlock(&hw_info_lock);
    if (ret != 0)
    {
        gmx_fatal(FARGS, "Error unlocking hwinfo mutex: %s", strerror(errno));
    }

    return hwinfo_g;
}

static std::string detected_hardware_string(const gmx_hw_info_t *hwinfo,
                                            bool                 bFullCpuInfo)
{
    std::string                  s;

    const gmx::CpuInfo          &cpuInfo = *hwinfo_g->cpuInfo;
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
    char host[HOSTNAMELEN];
    int  rank;

    gmx_gethostname(host, HOSTNAMELEN);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
        if (char *p = getenv("HWLOC_XMLFILE"))
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
            char buf[STRLEN];

            sprint_gpus(buf, &hwinfo->gpu_info);
            s += gmx::formatString("%s\n", buf);
        }
    }
    return s;
}

void gmx_print_detected_hardware(FILE *fplog, const t_commrec *cr,
                                 const gmx_hw_info_t *hwinfo)
{
    const gmx::CpuInfo &cpuInfo = *hwinfo_g->cpuInfo;

    if (fplog != NULL)
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
    check_use_of_rdtscp_on_this_cpu(fplog, cr, cpuInfo);
}

//! \brief Return if any GPU ID (e.g in a user-supplied string) is repeated
static gmx_bool anyGpuIdIsRepeated(const gmx_gpu_opt_t *gpu_opt)
{
    /* Loop over IDs in the string */
    for (int i = 0; i < gpu_opt->n_dev_use - 1; ++i)
    {
        /* Look for the ID in location i in the following part of the
           string */
        for (int j = i + 1; j < gpu_opt->n_dev_use; ++j)
        {
            if (gpu_opt->dev_use[i] == gpu_opt->dev_use[j])
            {
                /* Same ID found in locations i and j */
                return TRUE;
            }
        }
    }

    return FALSE;
}

void gmx_parse_gpu_ids(gmx_gpu_opt_t *gpu_opt)
{
    char *env;

    if (gpu_opt->gpu_id != NULL && !bGPUBinary)
    {
        gmx_fatal(FARGS, "GPU ID string set, but %s was compiled without GPU support!",
                  gmx::getProgramContext().displayName());
    }

    env = getenv("GMX_GPU_ID");
    if (env != NULL && gpu_opt->gpu_id != NULL)
    {
        gmx_fatal(FARGS, "GMX_GPU_ID and -gpu_id can not be used at the same time");
    }
    if (env == NULL)
    {
        env = gpu_opt->gpu_id;
    }

    /* parse GPU IDs if the user passed any */
    if (env != NULL)
    {
        /* Parse a "plain" or comma-separated GPU ID string which contains a
         * sequence of digits corresponding to GPU IDs; the order will
         * indicate the process/tMPI thread - GPU assignment. */
        parse_digits_from_string(env, &gpu_opt->n_dev_use, &gpu_opt->dev_use);

        if (!gmx_multiple_gpu_per_node_supported() && 1 < gpu_opt->n_dev_use)
        {
            gmx_fatal(FARGS, "The %s implementation only supports using exactly one PP rank per node", getGpuImplementationString());
        }
        if (!gmx_gpu_sharing_supported() && anyGpuIdIsRepeated(gpu_opt))
        {
            gmx_fatal(FARGS, "The %s implementation only supports using exactly one PP rank per GPU", getGpuImplementationString());
        }
        if (gpu_opt->n_dev_use == 0)
        {
            gmx_fatal(FARGS, "Empty GPU ID string encountered.\n%s\n",
                      invalid_gpuid_hint);
        }

        gpu_opt->bUserSet = TRUE;
    }
}

void gmx_select_gpu_ids(FILE *fplog, const t_commrec *cr,
                        const gmx_gpu_info_t *gpu_info,
                        gmx_bool bForceUseGPU,
                        gmx_gpu_opt_t *gpu_opt)
{
    int              i;
    char             sbuf[STRLEN], stmp[STRLEN];

    /* Bail if binary is not compiled with GPU acceleration, but this is either
     * explicitly (-nb gpu) or implicitly (gpu ID passed) requested. */
    if (bForceUseGPU && !bGPUBinary)
    {
        gmx_fatal(FARGS, "GPU acceleration requested, but %s was compiled without GPU support!",
                  gmx::getProgramContext().displayName());
    }

    if (!(cr->duty & DUTY_PP))
    {
        /* Our rank is not doing PP, we don't use a GPU */
        return;
    }

    if (gpu_opt->bUserSet)
    {
        /* Check the GPU IDs passed by the user.
         * (GPU IDs have been parsed by gmx_parse_gpu_ids before)
         */
        int *checkres;
        int  res;

        snew(checkres, gpu_opt->n_dev_use);

        res = check_selected_gpus(checkres, gpu_info, gpu_opt);

        if (!res)
        {
            print_gpu_detection_stats(fplog, gpu_info, cr);

            sprintf(sbuf, "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n");
            for (i = 0; i < gpu_opt->n_dev_use; i++)
            {
                if (checkres[i] != egpuCompatible)
                {
                    sprintf(stmp, "    GPU #%d: %s\n",
                            gpu_opt->dev_use[i],
                            gpu_detect_res_str[checkres[i]]);
                    strcat(sbuf, stmp);
                }
            }
            gmx_fatal(FARGS, "%s", sbuf);
        }

        sfree(checkres);
    }
    else if (getenv("GMX_EMULATE_GPU") == NULL)
    {
        pick_compatible_gpus(&hwinfo_g->gpu_info, gpu_opt);
        set_gpu_ids(gpu_opt, cr->nrank_pp_intranode, cr->rank_pp_intranode);
    }

    /* If the user asked for a GPU, check whether we have a GPU */
    if (bForceUseGPU && gpu_info->n_dev_compatible == 0)
    {
        gmx_fatal(FARGS, "GPU acceleration requested, but no compatible GPUs were detected.");
    }
}

/* Select the GPUs we will use. This is an operation local to each physical
 * node. If we have less MPI ranks than GPUs, we will waste some GPUs.
 * nrank and rank are the rank count and id for PP processes in our node.
 */
static void set_gpu_ids(gmx_gpu_opt_t *gpu_opt, int nrank, int rank)
{
    GMX_RELEASE_ASSERT(gpu_opt, "Invalid gpu_opt pointer passed");
    GMX_RELEASE_ASSERT(nrank >= 1,
                       gmx::formatString("Invalid limit (%d) for the number of GPUs (detected %d compatible GPUs)",
                                         rank, gpu_opt->n_dev_compatible).c_str());

    if (gpu_opt->n_dev_compatible == 0)
    {
        char host[HOSTNAMELEN];

        gmx_gethostname(host, HOSTNAMELEN);
        gmx_fatal(FARGS, "A GPU was requested on host %s, but no compatible GPUs were detected. All nodes with PP ranks need to have GPUs. If you intended to use GPU acceleration in a parallel run, you can either avoid using the nodes that don't have GPUs or place PME ranks on these nodes.", host);
    }

    int nshare;

    nshare = 1;
    if (nrank > gpu_opt->n_dev_compatible)
    {
        if (nrank % gpu_opt->n_dev_compatible == 0)
        {
            nshare = gmx_gpu_sharing_supported() ? nrank/gpu_opt->n_dev_compatible : 1;
        }
        else
        {
            if (rank == 0)
            {
                gmx_fatal(FARGS, "The number of MPI ranks (%d) in a physical node is not a multiple of the number of GPUs (%d). Select a different number of MPI ranks or use the -gpu_id option to manually specify the GPU to be used.",
                          nrank, gpu_opt->n_dev_compatible);
            }

#if GMX_MPI
            /* We use a global barrier to prevent ranks from continuing with
             * an invalid setup.
             */
            MPI_Barrier(MPI_COMM_WORLD);
#endif
        }
    }

    /* Here we will waste GPUs when nrank < gpu_opt->n_dev_compatible */
    gpu_opt->n_dev_use = std::min(gpu_opt->n_dev_compatible*nshare, nrank);
    if (!gmx_multiple_gpu_per_node_supported())
    {
        gpu_opt->n_dev_use = std::min(gpu_opt->n_dev_use, 1);
    }
    snew(gpu_opt->dev_use, gpu_opt->n_dev_use);
    for (int i = 0; i != gpu_opt->n_dev_use; ++i)
    {
        /* TODO: improve this implementation: either sort GPUs or remove the weakest here */
        gpu_opt->dev_use[i] = gpu_opt->dev_compatible[i/nshare];
    }
}

void gmx_hardware_info_free(gmx_hw_info_t *hwinfo)
{
    int ret;

    ret = tMPI_Thread_mutex_lock(&hw_info_lock);
    if (ret != 0)
    {
        gmx_fatal(FARGS, "Error locking hwinfo mutex: %s", strerror(errno));
    }

    /* decrease the reference counter */
    n_hwinfo--;


    if (hwinfo != hwinfo_g)
    {
        gmx_incons("hwinfo < hwinfo_g");
    }

    if (n_hwinfo < 0)
    {
        gmx_incons("n_hwinfo < 0");
    }

    if (n_hwinfo == 0)
    {
        delete hwinfo_g->cpuInfo;
        delete hwinfo_g->hardwareTopology;
        free_gpu_info(&hwinfo_g->gpu_info);
        sfree(hwinfo_g);
    }

    ret = tMPI_Thread_mutex_unlock(&hw_info_lock);
    if (ret != 0)
    {
        gmx_fatal(FARGS, "Error unlocking hwinfo mutex: %s", strerror(errno));
    }
}
