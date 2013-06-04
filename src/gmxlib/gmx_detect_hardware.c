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

/* Although we can't have more than 10 GPU different ID-s passed by the user as
 * the id-s are assumed to be represented by single digits, as multiple
 * processes can share a GPU, we can end up with more than 10 IDs.
 * To account for potential extreme cases we'll set the limit to a pretty
 * ridiculous number. */
static unsigned int max_gpu_ids_user = 64;

static const char * invalid_gpuid_hint =
    "A delimiter-free sequence of valid numeric IDs of available GPUs is expected.";

/* The globally shared hwinfo structure. */
static gmx_hw_info_t      *hwinfo_g;
/* A reference counter for the hwinfo structure */
static int                 n_hwinfo = 0;
/* A lock to protect the hwinfo structure */
static tMPI_Thread_mutex_t hw_info_lock = TMPI_THREAD_MUTEX_INITIALIZER;


/* FW decl. */
static void limit_num_gpus_used(gmx_hw_info_t *hwinfo, int count);

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
                                const t_commrec      *cr)
{
    char sbuf[STRLEN], stmp[STRLEN];
    int  i, ngpu, ngpu_all;

    ngpu     = gpu_info->ncuda_dev_use;
    ngpu_all = gpu_info->ncuda_dev;

    /* Issue note if GPUs are available but not used */
    if (ngpu_all > 0 && ngpu < 1)
    {
        sprintf(sbuf,
                "%d compatible GPU%s detected in the system, but none will be used.\n"
                "Consider trying GPU acceleration with the Verlet scheme!",
                ngpu_all, (ngpu_all > 1) ? "s" : "");
    }
    else
    {
        sprintf(sbuf, "%d GPU%s %sselected for this run: ",
                ngpu, (ngpu > 1) ? "s" : "",
                gpu_info->bUserSet ? "user-" : "auto-");
        for (i = 0; i < ngpu; i++)
        {
            sprintf(stmp, "#%d", get_gpu_device_id(gpu_info, i));
            if (i < ngpu - 1)
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
static void parse_gpu_id_plain_string(const char *idstr, int *nid, int *idlist)
{
    int    i;
    size_t len_idstr;

    len_idstr = strlen(idstr);

    if (len_idstr > max_gpu_ids_user)
    {
        gmx_fatal(FARGS, "%d GPU IDs provided, but only at most %d are supported",
                  len_idstr, max_gpu_ids_user);
    }

    *nid = len_idstr;

    for (i = 0; i < *nid; i++)
    {
        if (idstr[i] < '0' || idstr[i] > '9')
        {
            gmx_fatal(FARGS, "Invalid character in GPU ID string: '%c'\n%s\n",
                      idstr[i], invalid_gpuid_hint);
        }
        idlist[i] = idstr[i] - '0';
    }
}

static void parse_gpu_id_csv_string(const char *idstr, int *nid, int *idlist)
{
    /* XXX implement cvs format to support more than 10 different GPUs in a box. */
    gmx_incons("Not implemented yet");
}

void gmx_check_hw_runconf_consistency(FILE *fplog, gmx_hw_info_t *hwinfo,
                                      const t_commrec *cr, int ntmpi_requested,
                                      gmx_bool bUseGPU)
{
    int                        npppn, ntmpi_pp, ngpu;
    char                       sbuf[STRLEN], th_or_proc[STRLEN], th_or_proc_plural[STRLEN], pernode[STRLEN];
    char                       gpu_plural[2];
    gmx_bool                   bGPUBin, btMPI, bMPI, bMaxMpiThreadsSet, bNthreadsAuto, bEmulateGPU;
    int                        ret;
    static tMPI_Thread_mutex_t cons_lock = TMPI_THREAD_MUTEX_INITIALIZER;


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

    /* We run this function only once, but must make sure that all threads
       that are alive run this function, so they get consistent data. We
       achieve this by mutual exclusion and returning if the structure is
       already properly checked & set */
    ret = tMPI_Thread_mutex_lock(&cons_lock);
    if (ret != 0)
    {
        gmx_fatal(FARGS, "Error locking cons mutex: %s", strerror(errno));
    }

    if (!hwinfo->bConsistencyChecked)
    {
        btMPI         = bMPI = FALSE;
        bNthreadsAuto = FALSE;
#if defined(GMX_THREAD_MPI)
        btMPI         = TRUE;
        bNthreadsAuto = (ntmpi_requested < 1);
#elif defined(GMX_LIB_MPI)
        bMPI  = TRUE;
#endif

#ifdef GMX_GPU
        bGPUBin      = TRUE;
#else
        bGPUBin      = FALSE;
#endif

        /* GPU emulation detection is done later, but we need here as well
         * -- uncool, but there's no elegant workaround */
        bEmulateGPU       = (getenv("GMX_EMULATE_GPU") != NULL);
        bMaxMpiThreadsSet = (getenv("GMX_MAX_MPI_THREADS") != NULL);

        /* check the acceleration mdrun is compiled with against hardware
           capabilities */
        /* TODO: Here we assume homogeneous hardware which is not necessarily
                 the case! Might not hurt to add an extra check over MPI. */
        gmx_cpuid_acceleration_check(hwinfo->cpuid_info, fplog);

        /* Need to ensure that we have enough GPUs:
         * - need one GPU per PP node
         * - no GPU oversubscription with tMPI
         * => keep on the GPU support, otherwise turn off (or bail if forced)
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

        if (bGPUBin)
        {
            print_gpu_detection_stats(fplog, &hwinfo->gpu_info, cr);
        }

        if (bUseGPU && hwinfo->bCanUseGPU && !bEmulateGPU)
        {
            ngpu = hwinfo->gpu_info.ncuda_dev_use;
            sprintf(gpu_plural, "%s", (ngpu > 1) ? "s" : "");

            /* number of tMPI threads atuo-adjusted */
            if (btMPI && bNthreadsAuto)
            {
                if (npppn < ngpu)
                {
                    if (hwinfo->gpu_info.bUserSet)
                    {
                        /* The user manually provided more GPUs than threads we
                           could automatically start. */
                        gmx_fatal(FARGS,
                                  "%d GPU%s provided, but only %d PP thread-MPI thread%s coud be started.\n"
                                  "%s requires one PP tread-MPI thread per GPU; use fewer GPUs%s.",
                                  ngpu, gpu_plural, npppn, th_or_proc_plural,
                                  ShortProgram(), bMaxMpiThreadsSet ? "\nor allow more threads to be used" : "");
                    }
                    else
                    {
                        /* There are more GPUs than tMPI threads; we have to
                           limit the number GPUs used. */
                        md_print_warn(cr, fplog,
                                      "NOTE: %d GPU%s were detected, but only %d PP thread-MPI thread%s can be started.\n"
                                      "      %s can use one GPU per PP tread-MPI thread, so only %d GPU%s will be used.%s\n",
                                      ngpu, gpu_plural, npppn,
                                      th_or_proc_plural,
                                      ShortProgram(), npppn,
                                      npppn > 1 ? "s" : "",
                                      bMaxMpiThreadsSet ? "\n      Also, you can allow more threads to be used by increasing GMX_MAX_MPI_THREADS" : "");

                        if (cr->rank_pp_intranode == 0)
                        {
                            limit_num_gpus_used(hwinfo, npppn);
                            ngpu = hwinfo->gpu_info.ncuda_dev_use;
                            sprintf(gpu_plural, "%s", (ngpu > 1) ? "s" : "");
                        }
                    }
                }
            }

            if (ngpu != npppn)
            {
                if (hwinfo->gpu_info.bUserSet)
                {
                    gmx_fatal(FARGS,
                              "Incorrect launch configuration: mismatching number of PP %s%s and GPUs%s.\n"
                              "%s was started with %d PP %s%s%s, but you provided %d GPU%s.",
                              th_or_proc, btMPI ? "s" : "es", pernode,
                              ShortProgram(), npppn, th_or_proc,
                              th_or_proc_plural, pernode, ngpu, gpu_plural);
                }
                else
                {
                    if (ngpu > npppn)
                    {
                        md_print_warn(cr, fplog,
                                      "NOTE: potentially sub-optimal launch configuration, %s started with less\n"
                                      "      PP %s%s%s than GPU%s available.\n"
                                      "      Each PP %s can use only one GPU, %d GPU%s%s will be used.\n",
                                      ShortProgram(), th_or_proc,
                                      th_or_proc_plural, pernode, gpu_plural,
                                      th_or_proc, npppn, gpu_plural, pernode);

                        if (bMPI || (btMPI && cr->rank_pp_intranode == 0))
                        {
                            limit_num_gpus_used(hwinfo, npppn);
                            ngpu = hwinfo->gpu_info.ncuda_dev_use;
                            sprintf(gpu_plural, "%s", (ngpu > 1) ? "s" : "");
                        }
                    }
                    else
                    {
                        /* Avoid duplicate error messages.
                         * Unfortunately we can only do this at the physical node
                         * level, since the hardware setup and MPI process count
                         * might be differ over physical nodes.
                         */
                        if (cr->rank_pp_intranode == 0)
                        {
                            gmx_fatal(FARGS,
                                      "Incorrect launch configuration: mismatching number of PP %s%s and GPUs%s.\n"
                                      "%s was started with %d PP %s%s%s, but only %d GPU%s were detected.",
                                      th_or_proc, btMPI ? "s" : "es", pernode,
                                      ShortProgram(), npppn, th_or_proc,
                                      th_or_proc_plural, pernode, ngpu,
                                      gpu_plural);
                        }
                    }
                }
            }

            {
                int      same_count;

                same_count = gmx_count_gpu_dev_shared(&(hwinfo->gpu_info));

                if (btMPI && same_count > 0)
                {
                    gmx_fatal(FARGS,
                              "Invalid GPU assignment: can't share a GPU among multiple thread-MPI threads.\n"
                              "Use MPI if you are sure that you want to assign GPU to multiple threads.");
                }

                if (same_count > 0)
                {
                    md_print_warn(cr, fplog,
                                  "NOTE: Potentially sub-optimal launch configuration: you assigned %s to\n"
                                  "      multiple %s%s; this should be avoided as it can cause\n"
                                  "      performance loss.\n",
                                  same_count > 1 ? "GPUs" : "a GPU", th_or_proc, btMPI ? "s" : "es");
                }
            }
            print_gpu_use_stats(fplog, &hwinfo->gpu_info, cr);
        }
        hwinfo->bConsistencyChecked = TRUE;
    }

    ret = tMPI_Thread_mutex_unlock(&cons_lock);
    if (ret != 0)
    {
        gmx_fatal(FARGS, "Error unlocking cons mutex: %s", strerror(errno));
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

int gmx_count_gpu_dev_shared(const gmx_gpu_info_t *gpu_info)
{
    int      same_count    = 0;
    int      ngpu          = gpu_info->ncuda_dev_use;

    if (gpu_info->bUserSet)
    {
        int      i, j;

        for (i = 0; i < ngpu - 1; i++)
        {
            for (j = i + 1; j < ngpu; j++)
            {
                same_count      += (gpu_info->cuda_dev_use[i] ==
                                    gpu_info->cuda_dev_use[j]);
            }
        }
    }

    return same_count;
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

gmx_hw_info_t *gmx_detect_hardware(FILE *fplog, const t_commrec *cr,
                                   gmx_bool bForceUseGPU, gmx_bool bTryUseGPU,
                                   const char *gpu_id)
{
    int              i;
    const char      *env;
    char             sbuf[STRLEN], stmp[STRLEN];
    gmx_hw_info_t   *hw;
    gmx_gpu_info_t   gpuinfo_auto, gpuinfo_user;
    gmx_bool         bGPUBin;
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
        hwinfo_g->bConsistencyChecked = FALSE;

        /* detect CPUID info; no fuss, we don't detect system-wide
         * -- sloppy, but that's it for now */
        if (gmx_cpuid_init(&hwinfo_g->cpuid_info) != 0)
        {
            gmx_fatal_collective(FARGS, cr, NULL, "CPUID detection failed!");
        }

        /* detect number of hardware threads */
        hwinfo_g->nthreads_hw_avail = get_nthreads_hw_avail(fplog, cr);

        /* detect GPUs */
        hwinfo_g->gpu_info.ncuda_dev_use  = 0;
        hwinfo_g->gpu_info.cuda_dev_use   = NULL;
        hwinfo_g->gpu_info.ncuda_dev      = 0;
        hwinfo_g->gpu_info.cuda_dev       = NULL;

#ifdef GMX_GPU
        bGPUBin      = TRUE;
#else
        bGPUBin      = FALSE;
#endif

        /* Bail if binary is not compiled with GPU acceleration, but this is either
         * explicitly (-nb gpu) or implicitly (gpu ID passed) requested. */
        if (bForceUseGPU && !bGPUBin)
        {
            gmx_fatal(FARGS, "GPU acceleration requested, but %s was compiled without GPU support!", ShortProgram());
        }
        if (gpu_id != NULL && !bGPUBin)
        {
            gmx_fatal(FARGS, "GPU ID string set, but %s was compiled without GPU support!", ShortProgram());
        }

        /* run the detection if the binary was compiled with GPU support */
        if (bGPUBin && getenv("GMX_DISABLE_GPU_DETECTION") == NULL)
        {
            char detection_error[STRLEN];

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

        if (bForceUseGPU || bTryUseGPU)
        {
            env = getenv("GMX_GPU_ID");
            if (env != NULL && gpu_id != NULL)
            {
                gmx_fatal(FARGS, "GMX_GPU_ID and -gpu_id can not be used at the same time");
            }
            if (env == NULL)
            {
                env = gpu_id;
            }

            /* parse GPU IDs if the user passed any */
            if (env != NULL)
            {
                int *gpuid, *checkres;
                int  nid, res;

                snew(gpuid, max_gpu_ids_user);
                snew(checkres, max_gpu_ids_user);

                parse_gpu_id_plain_string(env, &nid, gpuid);

                if (nid == 0)
                {
                    gmx_fatal(FARGS, "Empty GPU ID string encountered.\n%s\n",
                              invalid_gpuid_hint);
                }

                res = check_select_cuda_gpus(checkres, &hwinfo_g->gpu_info,
                                             gpuid, nid);

                if (!res)
                {
                    print_gpu_detection_stats(fplog, &hwinfo_g->gpu_info, cr);

                    sprintf(sbuf, "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n");
                    for (i = 0; i < nid; i++)
                    {
                        if (checkres[i] != egpuCompatible)
                        {
                            sprintf(stmp, "    GPU #%d: %s\n",
                                    gpuid[i], gpu_detect_res_str[checkres[i]]);
                            strcat(sbuf, stmp);
                        }
                    }
                    gmx_fatal(FARGS, "%s", sbuf);
                }

                hwinfo_g->gpu_info.bUserSet = TRUE;

                sfree(gpuid);
                sfree(checkres);
            }
            else
            {
                pick_compatible_gpus(&hwinfo_g->gpu_info);
                hwinfo_g->gpu_info.bUserSet = FALSE;
            }

            /* decide whether we can use GPU */
            hwinfo_g->bCanUseGPU = (hwinfo_g->gpu_info.ncuda_dev_use > 0);
            if (!hwinfo_g->bCanUseGPU && bForceUseGPU)
            {
                gmx_fatal(FARGS, "GPU acceleration requested, but no compatible GPUs were detected.");
            }
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

static void limit_num_gpus_used(gmx_hw_info_t *hwinfo, int count)
{
    int ndev_use;

    assert(hwinfo);

    ndev_use = hwinfo->gpu_info.ncuda_dev_use;

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
    hwinfo->gpu_info.ncuda_dev_use = count;
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
