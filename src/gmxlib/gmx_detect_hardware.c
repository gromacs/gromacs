/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "types/enums.h"
#include "types/hw_info.h"
#include "types/commrec.h"
#include "gmx_fatal.h"
#include "gmx_fatal_collective.h"
#include "smalloc.h"
#include "gpu_utils.h"
#include "statutil.h"
#include "gmx_detect_hardware.h"
#include "main.h"
#include "md_logging.h"

#include "thread_mpi/threads.h"

#if ((defined(WIN32) || defined( _WIN32 ) || defined(WIN64) || defined( _WIN64 )) && !(defined (__CYGWIN__) || defined (__CYGWIN32__)))
#include "windows.h"
#endif

#ifdef GMX_GPU
const gmx_bool bGPUBinary = TRUE;
#else
const gmx_bool bGPUBinary = FALSE;
#endif

static const char * invalid_gpuid_hint =
    "A delimiter-free sequence of valid numeric IDs of available GPUs is expected.";

/* The globally shared hwinfo structure. */
static gmx_hw_info_t      *hwinfo_g;
/* A reference counter for the hwinfo structure */
static int                 n_hwinfo = 0;
/* A lock to protect the hwinfo structure */
static tMPI_Thread_mutex_t hw_info_lock = TMPI_THREAD_MUTEX_INITIALIZER;


/* FW decl. */
static void limit_num_gpus_used(gmx_gpu_opt_t *gpu_opt, int count);
static int gmx_count_gpu_dev_unique(const gmx_gpu_info_t *gpu_info,
                                    const gmx_gpu_opt_t  *gpu_opt);

static void sprint_gpus(char *sbuf, const gmx_gpu_info_t *gpu_info, gmx_bool bPrintAll)
{
    int      i, ndev;
    char     stmp[STRLEN];

    ndev = gpu_info->ncuda_dev;

    sbuf[0] = '\0';
    for (i = 0; i < ndev; i++)
    {
        get_gpu_device_info_string(stmp, gpu_info, i);
        strcat(sbuf, "  ");
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
    char onhost[266], stmp[STRLEN];
    int  ngpu;

    if (!gpu_info->bDetectGPUs)
    {
        /* We skipped the detection, so don't print detection stats */
        return;
    }

    ngpu = gpu_info->ncuda_dev;

#if defined GMX_MPI && !defined GMX_THREAD_MPI
    /* We only print the detection on one, of possibly multiple, nodes */
    strncpy(onhost, " on host ", 10);
    gmx_gethostname(onhost+9, 256);
#else
    /* We detect all relevant GPUs */
    strncpy(onhost, "", 1);
#endif

    if (ngpu > 0)
    {
        sprint_gpus(stmp, gpu_info, TRUE);
        md_print_warn(cr, fplog, "%d GPU%s detected%s:\n%s\n",
                      ngpu, (ngpu > 1) ? "s" : "", onhost, stmp);
    }
    else
    {
        md_print_warn(cr, fplog, "No GPUs detected%s\n", onhost);
    }
}

static void print_gpu_use_stats(FILE                 *fplog,
                                const gmx_gpu_info_t *gpu_info,
                                const gmx_gpu_opt_t  *gpu_opt,
                                const t_commrec      *cr)
{
    char sbuf[STRLEN], stmp[STRLEN];
    int  i, ngpu_comp, ngpu_use;

    ngpu_comp = gpu_info->ncuda_dev_compatible;
    ngpu_use  = gpu_opt->ncuda_dev_use;

    /* Issue a note if GPUs are available but not used */
    if (ngpu_comp > 0 && ngpu_use < 1)
    {
        sprintf(sbuf,
                "%d compatible GPU%s detected in the system, but none will be used.\n"
                "Consider trying GPU acceleration with the Verlet scheme!",
                ngpu_comp, (ngpu_comp > 1) ? "s" : "");
    }
    else
    {
        int ngpu_use_uniq;

        ngpu_use_uniq = gmx_count_gpu_dev_unique(gpu_info, gpu_opt);

        sprintf(sbuf, "%d GPU%s %sselected for this run.\n"
                "Mapping of GPU%s to the %d PP rank%s in this node: ",
                ngpu_use_uniq, (ngpu_use_uniq > 1) ? "s" : "",
                gpu_opt->bUserSet ? "user-" : "auto-",
                (ngpu_use > 1) ? "s" : "",
                cr->nrank_pp_intranode,
                (cr->nrank_pp_intranode > 1) ? "s" : "");

        for (i = 0; i < ngpu_use; i++)
        {
            sprintf(stmp, "#%d", get_gpu_device_id(gpu_info, gpu_opt, i));
            if (i < ngpu_use - 1)
            {
                strcat(stmp, ", ");
            }
            strcat(sbuf, stmp);
        }
    }
    md_print_info(cr, fplog, "%s\n\n", sbuf);
}

/* Parse a "plain" GPU ID string which contains a sequence of digits corresponding
 * to GPU IDs; the order will indicate the process/tMPI thread - GPU assignment. */
static void parse_gpu_id_plain_string(const char *idstr, int *nid, int **idlist)
{
    int i;

    *nid = strlen(idstr);

    snew(*idlist, *nid);

    for (i = 0; i < *nid; i++)
    {
        if (idstr[i] < '0' || idstr[i] > '9')
        {
            gmx_fatal(FARGS, "Invalid character in GPU ID string: '%c'\n%s\n",
                      idstr[i], invalid_gpuid_hint);
        }
        (*idlist)[i] = idstr[i] - '0';
    }
}

static void parse_gpu_id_csv_string(const char *idstr, int *nid, int *idlist)
{
    /* XXX implement cvs format to support more than 10 different GPUs in a box. */
    gmx_incons("Not implemented yet");
}

void gmx_check_hw_runconf_consistency(FILE *fplog,
                                      const gmx_hw_info_t *hwinfo,
                                      const t_commrec *cr,
                                      const gmx_hw_opt_t *hw_opt,
                                      gmx_bool bUseGPU)
{
    int      npppn, ntmpi_pp;
    char     sbuf[STRLEN], th_or_proc[STRLEN], th_or_proc_plural[STRLEN], pernode[STRLEN];
    gmx_bool btMPI, bMPI, bMaxMpiThreadsSet, bNthreadsAuto, bEmulateGPU;

    assert(hwinfo);
    assert(cr);

    /* Below we only do consistency checks for PP and GPUs,
     * this is irrelevant for PME only nodes, so in that case we return
     * here.
     */
    if (!(cr->duty & DUTY_PP))
    {
        return;
    }

    btMPI         = bMPI = FALSE;
    bNthreadsAuto = FALSE;
#if defined(GMX_THREAD_MPI)
    btMPI         = TRUE;
    bNthreadsAuto = (hw_opt->nthreads_tmpi < 1);
#elif defined(GMX_LIB_MPI)
    bMPI  = TRUE;
#endif

    /* GPU emulation detection is done later, but we need here as well
     * -- uncool, but there's no elegant workaround */
    bEmulateGPU       = (getenv("GMX_EMULATE_GPU") != NULL);
    bMaxMpiThreadsSet = (getenv("GMX_MAX_MPI_THREADS") != NULL);

    /* check the acceleration mdrun is compiled with against hardware
       capabilities */
    /* TODO: Here we assume homogeneous hardware which is not necessarily
             the case! Might not hurt to add an extra check over MPI. */
    gmx_cpuid_acceleration_check(hwinfo->cpuid_info, fplog, SIMMASTER(cr));

    /* NOTE: this print is only for and on one physical node */
    print_gpu_detection_stats(fplog, &hwinfo->gpu_info, cr);

    if (hwinfo->gpu_info.ncuda_dev_compatible > 0)
    {
        /* NOTE: this print is only for and on one physical node */
        print_gpu_use_stats(fplog, &hwinfo->gpu_info, &hw_opt->gpu_opt, cr);
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

    if (bUseGPU && hwinfo->gpu_info.ncuda_dev_compatible > 0 &&
        !bEmulateGPU)
    {
        int  ngpu_comp, ngpu_use;
        char gpu_comp_plural[2], gpu_use_plural[2];

        ngpu_comp = hwinfo->gpu_info.ncuda_dev_compatible;
        ngpu_use  = hw_opt->gpu_opt.ncuda_dev_use;

        sprintf(gpu_comp_plural, "%s", (ngpu_comp> 1) ? "s" : "");
        sprintf(gpu_use_plural,  "%s", (ngpu_use > 1) ? "s" : "");

        /* number of tMPI threads auto-adjusted */
        if (btMPI && bNthreadsAuto)
        {
            if (hw_opt->gpu_opt.bUserSet && npppn < ngpu_use)
            {
                /* The user manually provided more GPUs than threads we
                   could automatically start. */
                gmx_fatal(FARGS,
                          "%d GPU%s provided, but only %d PP thread-MPI thread%s coud be started.\n"
                          "%s requires one PP tread-MPI thread per GPU; use fewer GPUs%s.",
                          ngpu_use, gpu_use_plural,
                          npppn, th_or_proc_plural,
                          ShortProgram(), bMaxMpiThreadsSet ? "\nor allow more threads to be used" : "");
            }

            if (!hw_opt->gpu_opt.bUserSet && npppn < ngpu_comp)
            {
                /* There are more GPUs than tMPI threads; we have
                   limited the number GPUs used. */
                md_print_warn(cr, fplog,
                              "NOTE: %d GPU%s were detected, but only %d PP thread-MPI thread%s can be started.\n"
                              "      %s can use one GPU per PP tread-MPI thread, so only %d GPU%s will be used.%s\n",
                              ngpu_comp, gpu_comp_plural,
                              npppn, th_or_proc_plural,
                              ShortProgram(), npppn,
                              npppn > 1 ? "s" : "",
                              bMaxMpiThreadsSet ? "\n      Also, you can allow more threads to be used by increasing GMX_MAX_MPI_THREADS" : "");
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
                          ShortProgram(), npppn, th_or_proc,
                          th_or_proc_plural, pernode,
                          ngpu_use, gpu_use_plural);
            }
        }
        else
        {
            if (ngpu_comp > npppn)
            {
                md_print_warn(cr, fplog,
                              "NOTE: potentially sub-optimal launch configuration, %s started with less\n"
                              "      PP %s%s%s than GPU%s available.\n"
                              "      Each PP %s can use only one GPU, %d GPU%s%s will be used.\n",
                              ShortProgram(), th_or_proc,
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
                    gmx_fatal(FARGS,
                              "Incorrect launch configuration: mismatching number of PP %s%s and GPUs%s.\n"
                              "%s was started with %d PP %s%s%s, but only %d GPU%s were detected.",
                              th_or_proc, btMPI ? "s" : "es", pernode,
                              ShortProgram(), npppn, th_or_proc,
                              th_or_proc_plural, pernode,
                              ngpu_use, gpu_use_plural);
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

#ifdef GMX_MPI
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
    int      ngpu          = gpu_opt->ncuda_dev_use;

    if (gpu_opt->bUserSet)
    {
        int      i, j;

        for (i = 0; i < ngpu - 1; i++)
        {
            for (j = i + 1; j < ngpu; j++)
            {
                same_count      += (gpu_opt->cuda_dev_use[i] ==
                                    gpu_opt->cuda_dev_use[j]);
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

    assert(gpu_info);
    assert(gpu_opt);

    ngpu        = gpu_info->ncuda_dev;
    uniq_count  = 0;

    snew(uniq_ids, ngpu);

    /* Each element in uniq_ids will be set to 0 or 1. The n-th element set
        * to 1 indicates that the respective GPU was selected to be used. */
    for (i = 0; i < gpu_opt->ncuda_dev_use; i++)
    {
        uniq_ids[get_gpu_device_id(gpu_info, gpu_opt, i)] = 1;
    }
    /* Count the devices used. */
    for (i = 0; i < ngpu; i++)
    {
        uniq_count += uniq_ids[i];
    }

    sfree(uniq_ids);

    return uniq_count;
}


/* Return the number of hardware threads supported by the current CPU.
 * We assume that this is equal with the number of CPUs reported to be
 * online by the OS at the time of the call.
 */
static int get_nthreads_hw_avail(FILE *fplog, const t_commrec *cr)
{
    int ret = 0;

#if ((defined(WIN32) || defined( _WIN32 ) || defined(WIN64) || defined( _WIN64 )) && !(defined (__CYGWIN__) || defined (__CYGWIN32__)))
    /* Windows */
    SYSTEM_INFO sysinfo;
    GetSystemInfo( &sysinfo );
    ret = sysinfo.dwNumberOfProcessors;
#elif defined HAVE_SYSCONF
    /* We are probably on Unix.
     * Now check if we have the argument to use before executing the call
     */
#if defined(_SC_NPROCESSORS_ONLN)
    ret = sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROC_ONLN)
    ret = sysconf(_SC_NPROC_ONLN);
#elif defined(_SC_NPROCESSORS_CONF)
    ret = sysconf(_SC_NPROCESSORS_CONF);
#elif defined(_SC_NPROC_CONF)
    ret = sysconf(_SC_NPROC_CONF);
#endif /* End of check for sysconf argument values */

#else
    /* Neither windows nor Unix. No fscking idea how many CPUs we have! */
    ret = -1;
#endif

    if (debug)
    {
        fprintf(debug, "Detected %d processors, will use this as the number "
                "of supported hardware threads.\n", ret);
    }

#ifdef GMX_OMPENMP
    if (ret != gmx_omp_get_num_procs())
    {
        md_print_warn(cr, fplog,
                      "Number of CPUs detected (%d) does not match the number reported by OpenMP (%d).\n"
                      "Consider setting the launch configuration manually!",
                      ret, gmx_omp_get_num_procs());
    }
#endif

    return ret;
}

static void gmx_detect_gpus(FILE *fplog, const t_commrec *cr,
                            gmx_gpu_info_t *gpu_info)
{
#ifdef GMX_LIB_MPI
    int              rank_world;
    MPI_Comm         physicalnode_comm;
#endif
    int              rank_local;

    /* Under certain circumstances MPI ranks on the same physical node
     * can not simultaneously access the same GPU(s). Therefore we run
     * the detection only on one MPI rank per node and broadcast the info.
     * Note that with thread-MPI only a single thread runs this code.
     *
     * TODO: We should also do CPU hardware detection only once on each
     * physical node and broadcast it, instead of do it on every MPI rank.
     */
#ifdef GMX_LIB_MPI
    /* A split of MPI_COMM_WORLD over physical nodes is only required here,
     * so we create and destroy it locally.
     */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
    MPI_Comm_split(MPI_COMM_WORLD, gmx_physicalnode_id_hash(),
                   rank_world, &physicalnode_comm);
    MPI_Comm_rank(physicalnode_comm, &rank_local);
#else
    /* Here there should be only one process, check this */
    assert(cr->nnodes == 1 && cr->sim_nodeid == 0);

    rank_local = 0;
#endif

    if (rank_local == 0)
    {
        char detection_error[STRLEN], sbuf[STRLEN];

        if (detect_cuda_gpus(&hwinfo_g->gpu_info, detection_error) != 0)
        {
            if (detection_error != NULL && detection_error[0] != '\0')
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

#ifdef GMX_LIB_MPI
    /* Broadcast the GPU info to the other ranks within this node */
    MPI_Bcast(&hwinfo_g->gpu_info.ncuda_dev, 1, MPI_INT, 0, physicalnode_comm);

    if (hwinfo_g->gpu_info.ncuda_dev > 0)
    {
        int cuda_dev_size;

        cuda_dev_size = hwinfo_g->gpu_info.ncuda_dev*sizeof_cuda_dev_info();

        if (rank_local > 0)
        {
            hwinfo_g->gpu_info.cuda_dev =
                (cuda_dev_info_ptr_t)malloc(cuda_dev_size);
        }
        MPI_Bcast(hwinfo_g->gpu_info.cuda_dev, cuda_dev_size, MPI_BYTE,
                  0, physicalnode_comm);
        MPI_Bcast(&hwinfo_g->gpu_info.ncuda_dev_compatible, 1, MPI_INT,
                  0, physicalnode_comm);
    }

    MPI_Comm_free(&physicalnode_comm);
#endif
}

gmx_hw_info_t *gmx_detect_hardware(FILE *fplog, const t_commrec *cr,
                                   gmx_bool bDetectGPUs)
{
    gmx_hw_info_t   *hw;
    int              ret;

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

        /* detect CPUID info; no fuss, we don't detect system-wide
         * -- sloppy, but that's it for now */
        if (gmx_cpuid_init(&hwinfo_g->cpuid_info) != 0)
        {
            gmx_fatal_collective(FARGS, cr, NULL, "CPUID detection failed!");
        }

        /* detect number of hardware threads */
        hwinfo_g->nthreads_hw_avail = get_nthreads_hw_avail(fplog, cr);

        /* detect GPUs */
        hwinfo_g->gpu_info.ncuda_dev            = 0;
        hwinfo_g->gpu_info.cuda_dev             = NULL;
        hwinfo_g->gpu_info.ncuda_dev_compatible = 0;

        /* Run the detection if the binary was compiled with GPU support
         * and we requested detection.
         */
        hwinfo_g->gpu_info.bDetectGPUs =
            (bGPUBinary && bDetectGPUs &&
             getenv("GMX_DISABLE_GPU_DETECTION") == NULL);
        if (hwinfo_g->gpu_info.bDetectGPUs)
        {
            gmx_detect_gpus(fplog, cr, &hwinfo_g->gpu_info);
        }
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

void gmx_parse_gpu_ids(gmx_gpu_opt_t *gpu_opt)
{
    char *env;

    if (gpu_opt->gpu_id != NULL && !bGPUBinary)
    {
        gmx_fatal(FARGS, "GPU ID string set, but %s was compiled without GPU support!", ShortProgram());
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
        parse_gpu_id_plain_string(env,
                                  &gpu_opt->ncuda_dev_use,
                                  &gpu_opt->cuda_dev_use);

        if (gpu_opt->ncuda_dev_use == 0)
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
    const char      *env;
    char             sbuf[STRLEN], stmp[STRLEN];

    /* Bail if binary is not compiled with GPU acceleration, but this is either
     * explicitly (-nb gpu) or implicitly (gpu ID passed) requested. */
    if (bForceUseGPU && !bGPUBinary)
    {
        gmx_fatal(FARGS, "GPU acceleration requested, but %s was compiled without GPU support!", ShortProgram());
    }

    if (gpu_opt->bUserSet)
    {
        /* Check the GPU IDs passed by the user.
         * (GPU IDs have been parsed by gmx_parse_gpu_ids before)
         */
        int *checkres;
        int  res;

        snew(checkres, gpu_opt->ncuda_dev_use);

        res = check_selected_cuda_gpus(checkres, gpu_info, gpu_opt);

        if (!res)
        {
            print_gpu_detection_stats(fplog, gpu_info, cr);

            sprintf(sbuf, "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n");
            for (i = 0; i < gpu_opt->ncuda_dev_use; i++)
            {
                if (checkres[i] != egpuCompatible)
                {
                    sprintf(stmp, "    GPU #%d: %s\n",
                            gpu_opt->cuda_dev_use[i],
                            gpu_detect_res_str[checkres[i]]);
                    strcat(sbuf, stmp);
                }
            }
            gmx_fatal(FARGS, "%s", sbuf);
        }

        sfree(checkres);
    }
    else
    {
        pick_compatible_gpus(&hwinfo_g->gpu_info, gpu_opt);

        if (gpu_opt->ncuda_dev_use > cr->nrank_pp_intranode)
        {
            /* We picked more GPUs than we can use: limit the number.
             * We print detailed messages about this later in
             * gmx_check_hw_runconf_consistency.
             */
            limit_num_gpus_used(gpu_opt, cr->nrank_pp_intranode);
        }

        gpu_opt->bUserSet = FALSE;
    }

    /* If the user asked for a GPU, check whether we have a GPU */
    if (bForceUseGPU && gpu_info->ncuda_dev_compatible == 0)
    {
        gmx_fatal(FARGS, "GPU acceleration requested, but no compatible GPUs were detected.");
    }
}

static void limit_num_gpus_used(gmx_gpu_opt_t *gpu_opt, int count)
{
    int ndev_use;

    assert(gpu_opt);

    ndev_use = gpu_opt->ncuda_dev_use;

    if (count > ndev_use)
    {
        /* won't increase the # of GPUs */
        return;
    }

    if (count < 1)
    {
        char sbuf[STRLEN];
        sprintf(sbuf, "Limiting the number of GPUs to <1 doesn't make sense (detected %d, %d requested)!",
                ndev_use, count);
        gmx_incons(sbuf);
    }

    /* TODO: improve this implementation: either sort GPUs or remove the weakest here */
    gpu_opt->ncuda_dev_use = count;
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
        gmx_cpuid_done(hwinfo_g->cpuid_info);
        free_gpu_info(&hwinfo_g->gpu_info);
        sfree(hwinfo_g);
    }

    ret = tMPI_Thread_mutex_unlock(&hw_info_lock);
    if (ret != 0)
    {
        gmx_fatal(FARGS, "Error unlocking hwinfo mutex: %s", strerror(errno));
    }
}
