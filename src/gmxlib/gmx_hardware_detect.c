/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of GROMACS.
 * Copyright (c) 2012-  
 *
 * Written by the Gromacs development team under coordination of
 * David van der Spoel, Berk Hess, and Erik Lindahl.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * GROup of MAchos and Cynical Suckers
 */

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "types/enums.h"
#include "types/hwinfo.h"
#include "types/commrec.h"
#include "gmx_fatal.h"
#include "gmx_fatal_collective.h"
#include "smalloc.h"
#include "gpu_utils.h"
#include "statutil.h"
#include "gmx_hardware_detect.h"
#include "main.h"
#include "md_logging.h"


/* FW decl. */
void limit_num_gpus_used(gmx_hwinfo_t *hwinfo, int count);

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

static void print_gpu_detection_stats(FILE *fplog,
                                      const gmx_gpu_info_t *gpu_info,
                                      const t_commrec *cr)
{
    char onhost[266],stmp[STRLEN];
    int  ngpu;

    ngpu = gpu_info->ncuda_dev;

#if defined GMX_MPI && !defined GMX_THREAD_MPI
    /* We only print the detection on one, of possibly multiple, nodes */
    strncpy(onhost," on host ",10);
    gmx_gethostname(onhost+9,256);
#else
    /* We detect all relevant GPUs */
    strncpy(onhost,"",1);
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

static void print_gpu_use_stats(FILE *fplog,
                                const gmx_gpu_info_t *gpu_info,
                                const t_commrec *cr)
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
        sprintf(sbuf, "%d GPU%s %sselected to be used for this run: ",
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

/* We can't have more than 10 GPU id's as we use single char numbers */
#define MAX_GPU_IDS  10

/* Parse a "plain" GPU ID string which contains a sequence of digits corresponding
 * to GPU IDs; the order will indicate the process/tMPI thread - GPU assignment. */
static void parse_gpu_id_plain_string(const char *idstr, int *nid, int *idlist)
{
    int  i;

    *nid = strlen(idstr);

    if (*nid > MAX_GPU_IDS)
    {
        gmx_fatal(FARGS,"%d GPU id's passed in a string, but not more than %d are supported",*nid,MAX_GPU_IDS);
    }

    for (i = 0; i < *nid; i++)
    {
        if (idstr[i] < '0' || idstr[i] > '9')
        {
            gmx_fatal(FARGS, "Invalid character in GPU ID string: '%c'\n", idstr[i]);
        }
        idlist[i] = idstr[i] - '0';
    }
}

static void parse_gpu_id_csv_string(const char *idstr, int *nid, int *idlist)
{
    /* TODO implement */
    gmx_incons("Not implemented yet");
}

void gmx_check_hw_runconf_consistency(FILE *fplog, gmx_hwinfo_t *hwinfo,
                                      const t_commrec *cr, int ntmpi_requested,
                                      gmx_bool bUseGPU)
{
    int      npppn, ntmpi_pp, ngpu;
    char     sbuf[STRLEN], th_or_proc[STRLEN], th_or_proc_plural[STRLEN], pernode[STRLEN];
    char     gpu_plural[2];
    gmx_bool bGPUBin, btMPI, bMPI, bMaxMpiThreadsSet, bNthreadsAuto, bEmulateGPU;

    assert(hwinfo);
    assert(cr);

    btMPI = bMPI = FALSE;
    bNthreadsAuto = FALSE;
#if defined(GMX_THREAD_MPI)
    btMPI = TRUE;
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

    if(SIMMASTER(cr))
    {
        /* check the acceleration mdrun is compiled with against hardware capabilities */
        /* TODO: Here we assume homogenous hardware which is not necessarily the case!
         *       Might not hurt to add an extra check over MPI. */
        gmx_detectcpu_check_acceleration(hwinfo->cpu_info, fplog);
    }

    /* Need to ensure that we have enough GPUs:
     * - need one GPU per PP node
     * - no GPU oversubscription with tMPI
     * => keep on the GPU support, otherwise turn off (or bail if forced)
     * */
    /* number of PP processes per node */
    npppn = cr->nnodes_pp_intra;

    pernode[0] = '\0';
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
        if (btMPI && bNthreadsAuto && SIMMASTER(cr))
        {
            if (npppn < ngpu)
            {
                if (hwinfo->gpu_info.bUserSet)
                {
                    /* The user manually provided more GPUs than threads we could
                     * automatically start. */
                    gmx_fatal(FARGS,
                              "%d GPU%s provided, but only %d PP thread-MPI thread%s coud be started.\n"
                              "%s requires one PP tread-MPI thread per GPU; use fewer GPUs%s.",
                              ngpu, gpu_plural, npppn, th_or_proc_plural,
                              ShortProgram(), bMaxMpiThreadsSet ? "\nor allow more threads to be used" : "");
                }
                else
                {
                    /* There are more GPUs than tMPI threads; we have to limit the number GPUs used. */
                    md_print_warn(cr,fplog,
                                  "NOTE: %d GPU%s were detected, but only %d PP thread-MPI thread%s can be started.\n"
                                  "      %s can use one GPU per PP tread-MPI thread, so only %d GPU%s will be used.%s\n",
                                  ngpu, gpu_plural, npppn, th_or_proc_plural,
                                  ShortProgram(), npppn, npppn > 1 ? "s" : "",
                                  bMaxMpiThreadsSet ? "\n      Also, you can allow more threads to be used by increasing GMX_MAX_MPI_THREADS" : "");

                    if (cr->nodeid_intra == 0)
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
                          th_or_proc, btMPI ? "s" : "es" , pernode,
                          ShortProgram(), npppn, th_or_proc, th_or_proc_plural, pernode, ngpu, gpu_plural);
            }
            else
            {
                if (ngpu > npppn)
                {
                    md_print_warn(cr,fplog,
                                  "NOTE: potentially sub-optimal launch configuration, %s started with less\n"
                                  "      PP %s%s%s than GPU%s available.\n"
                                  "      Each PP %s can only use one GPU, so only %d GPU%s%s will be used.",
                                  ShortProgram(),
                                  th_or_proc, th_or_proc_plural, pernode, gpu_plural,
                                  th_or_proc, npppn, gpu_plural, pernode);

                    if (bMPI || (btMPI && cr->nodeid_intra == 0))
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
                    if (cr->nodeid_intra == 0)
                    {
                        gmx_fatal(FARGS,
                                  "Incorrect launch configuration: mismatching number of PP %s%s and GPUs%s.\n"
                                  "%s was started with %d PP %s%s%s, but only %d GPU%s were detected.",
                                  th_or_proc, btMPI ? "s" : "es" , pernode,
                                  ShortProgram(), npppn, th_or_proc, th_or_proc_plural, pernode, ngpu, gpu_plural);
                    }
#ifdef GMX_MPI
                    else
                    {
                        /* Avoid other ranks to continue after inconsistency */
                        MPI_Barrier(MPI_COMM_WORLD);
                    }
#endif
                }
            }
        }

        if (hwinfo->gpu_info.bUserSet && (cr->nodeid_intra == 0))
        {
            int i, j, same_count;
            gmx_bool bSomeSame, bAllDifferent;

            same_count = 0;
            bSomeSame = FALSE;
            bAllDifferent = TRUE;

            for (i = 0; i < ngpu - 1; i++)
            {
                for (j = i + 1; j < ngpu; j++)
                {
                    bSomeSame       |= hwinfo->gpu_info.cuda_dev_use[i] == hwinfo->gpu_info.cuda_dev_use[j];
                    bAllDifferent   &= hwinfo->gpu_info.cuda_dev_use[i] != hwinfo->gpu_info.cuda_dev_use[j];
                    same_count      += hwinfo->gpu_info.cuda_dev_use[i] == hwinfo->gpu_info.cuda_dev_use[j];
                }
            }

            if (btMPI && !bAllDifferent)
            {
                gmx_fatal(FARGS,
                          "Invalid GPU assignment: can't share a GPU among multiple thread-MPI threads.\n"
                          "Use MPI if you are sure that you want to assign GPU to multiple threads.");
            }

            if (bSomeSame)
            {
                md_print_warn(cr,fplog,
                              "NOTE: Potentially sub-optimal launch configuration: you assigned %s to\n"
                              "      multiple %s%s; this should be avoided as it generally\n"
                              "      causes performance loss.",
                              same_count > 1 ? "GPUs" : "a GPU", th_or_proc, btMPI ? "s" : "es");
            }
        }
        print_gpu_use_stats(fplog, &hwinfo->gpu_info, cr);
    }
}

void gmx_hw_detect(FILE *fplog, gmx_hwinfo_t *hwinfo,
                   const t_commrec *cr,
                   gmx_bool bForceUseGPU, gmx_bool bTryUseGPU,
                   const char *gpu_id)
{
    int             i;
    const char      *env;
    char            sbuf[STRLEN], stmp[STRLEN];
    gmx_hwinfo_t    *hw;
    gmx_gpu_info_t  gpuinfo_auto, gpuinfo_user;
    gmx_bool        bGPUBin;

    assert(hwinfo);

    /* detect CPU; no fuss, we don't detect system-wide -- sloppy, but that's it for now */
    gmx_detectcpu(&hwinfo->cpu_info);

    /* detect GPUs */
    hwinfo->gpu_info.ncuda_dev_use  = 0;
    hwinfo->gpu_info.cuda_dev_use   = NULL;
    hwinfo->gpu_info.ncuda_dev      = 0;
    hwinfo->gpu_info.cuda_dev       = NULL;

#ifdef GMX_GPU
    bGPUBin      = TRUE;
#else
    bGPUBin      = FALSE;
#endif

    /* Bail if binary is not compiled with GPU on */
    if (bForceUseGPU && !bGPUBin)
    {
        gmx_fatal_collective(FARGS, cr, NULL, "GPU acceleration requested, but %s was compiled without GPU support!", ShortProgram());
    }

    /* run the detection if the binary was compiled with GPU support */
    if (bGPUBin && getenv("GMX_DISABLE_GPU_DETECTION")==NULL)
    {
        detect_cuda_gpus(&hwinfo->gpu_info);
    }

    if (bForceUseGPU || bTryUseGPU)
    {
        env = getenv("GMX_GPU_ID");
        if (env != NULL && gpu_id != NULL)
        {
            gmx_fatal(FARGS,"GMX_GPU_ID and -gpu_id can not be used at the same time");
        }
        if (env == NULL)
        {
            env = gpu_id;
        }

        /* parse GPU IDs if the user passed any */
        if (env != NULL)
        {
            int gpuid[MAX_GPU_IDS];
            int nid, res;
            int checkres[MAX_GPU_IDS];

            parse_gpu_id_plain_string(env, &nid, gpuid);

            if (nid == 0)
            {
                gmx_fatal(FARGS, "Empty GPU ID string passed\n");
            }

            res = check_select_cuda_gpus(checkres, &hwinfo->gpu_info, gpuid, nid);

            if (!res)
            {
                print_gpu_detection_stats(fplog, &hwinfo->gpu_info, cr);

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

            hwinfo->gpu_info.bUserSet = TRUE;
        }
        else
        {
            pick_compatible_gpus(&hwinfo->gpu_info);
            hwinfo->gpu_info.bUserSet = FALSE;
        }

        /* decide whether we can use GPU */
        hwinfo->bCanUseGPU = (hwinfo->gpu_info.ncuda_dev_use > 0);
        if (!hwinfo->bCanUseGPU && bForceUseGPU)
        {
            gmx_fatal(FARGS, "GPU acceleration requested, but no compatible GPUs were detected.");
        }
    }
}

void limit_num_gpus_used(gmx_hwinfo_t *hwinfo, int count)
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

void gmx_hw_info_free(gmx_hwinfo_t *hwinfo)
{
    if (hwinfo)
    {
        free_gpu_info(&hwinfo->gpu_info);
        sfree(hwinfo);
    }
}
