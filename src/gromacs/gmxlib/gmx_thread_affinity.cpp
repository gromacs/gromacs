/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/gmx_thread_affinity.h"

#include "config.h"

#include <tbb/tbb.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>

#ifdef HAVE_SCHED_AFFINITY
#  include <sched.h>
#  include <sys/syscall.h>
#endif

#include "thread_mpi/threads.h"

#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/gmx_cpuid.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/md_logging.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/types/hw_info.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

static int
get_thread_affinity_layout(FILE *fplog,
                           const t_commrec *cr,
                           const gmx_hw_info_t * hwinfo,
                           int nthreads,
                           int pin_offset, int * pin_stride,
                           const int **locality_order)
{
    int         nhwthreads, npkg, ncores, nhwthreads_per_core, rc;
    const int * pkg_id;
    const int * core_id;
    const int * hwthread_id;
    gmx_bool    bPickPinStride;

    if (pin_offset < 0)
    {
        gmx_fatal(FARGS, "Negative thread pinning offset requested");
    }
    if (*pin_stride < 0)
    {
        gmx_fatal(FARGS, "Negative thread pinning stride requested");
    }

    rc = gmx_cpuid_topology(hwinfo->cpuid_info, &nhwthreads, &npkg, &ncores,
                            &nhwthreads_per_core,
                            &pkg_id, &core_id, &hwthread_id, locality_order);

    if (rc != 0)
    {
        /* topology information not available or invalid, ignore it */
        nhwthreads      = hwinfo->nthreads_hw_avail;
        *locality_order = NULL;

        if (nhwthreads <= 0)
        {
            /* We don't know anything about the hardware, don't pin */
            md_print_warn(cr, fplog,
                          "NOTE: We don't know how many logical cores we have, will not pin threads");

            return -1;
        }
    }

    if (nthreads > nhwthreads)
    {
        /* We are oversubscribing, don't pin */
        md_print_warn(NULL, fplog,
                      "WARNING: Oversubscribing the CPU, will not pin threads");

        return -1;
    }

    if (pin_offset + nthreads > nhwthreads)
    {
        /* We are oversubscribing, don't pin */
        md_print_warn(NULL, fplog,
                      "WARNING: The requested pin offset is too large for the available logical cores,\n"
                      "         will not pin threads");

        return -1;
    }


    /* do we need to choose the pinning stride? */
    bPickPinStride = (*pin_stride == 0);

    if (bPickPinStride)
    {
        if (rc == 0 && pin_offset + nthreads*nhwthreads_per_core <= nhwthreads)
        {
            /* Put one thread on each physical core */
            *pin_stride = nhwthreads_per_core;
        }
        else
        {
            /* We don't know if we have SMT, and if we do, we don't know
             * if hw threads in the same physical core are consecutive.
             * Without SMT the pinning layout should not matter too much.
             * so we assume a consecutive layout and maximally spread out"
             * the threads at equal threads per core.
             * Note that IBM is the major non-x86 case with cpuid support
             * and probably threads are already pinned by the queuing system,
             * so we wouldn't end up here in the first place.
             */
            *pin_stride = (nhwthreads - pin_offset)/nthreads;
        }
    }
    else
    {
        /* Check the placement of the thread with the largest index to make sure
         * that the offset & stride doesn't cause pinning beyond the last hardware thread. */
        if (pin_offset + (nthreads-1)*(*pin_stride) >= nhwthreads)
        {
            /* We are oversubscribing, don't pin */
            md_print_warn(NULL, fplog,
                          "WARNING: The requested pinning stride is too large for the available logical cores,\n"
                          "         will not pin threads");

            return -1;
        }
    }

    if (fplog != NULL)
    {
        fprintf(fplog, "Pinning threads with a%s logical core stride of %d\n",
                bPickPinStride ? "n auto-selected" : " user-specified",
                *pin_stride);
    }

    return 0;
}

/* Set CPU affinity for current thread only. Can be important for performance.
   On some systems (e.g. Cray) CPU Affinity is set by default.
   But default assigning doesn't work (well) with only some ranks
   having threads. This causes very low performance.
   External tools have cumbersome syntax for setting affinity
   in the case that only some ranks have threads.
   Thus it is important that GROMACS sets the affinity internally
   if only PME is using threads.
 */
static void
gmx_set_thread_affinity(FILE                *fplog,
                        const t_commrec     *cr,
                        gmx_hw_opt_t        *hw_opt,
                        const gmx_hw_info_t *hwinfo)
{
    int        nth_affinity_set, thread0_id_node,
               nthread_local, nthread_node, nthread_hw_max, nphyscore;
    int        offset;
    const int *locality_order;
    int        rc;

    if (hw_opt->thread_affinity == threadaffOFF)
    {
        /* Nothing to do */
        return;
    }

    /* If the tMPI thread affinity setting is not supported encourage the user
     * to report it as it's either a bug or an exotic platform which we might
     * want to support. */
    if (tMPI_Thread_setaffinity_support() != TMPI_SETAFFINITY_SUPPORT_YES)
    {
        /* we know Mac OS doesn't support setting thread affinity, so there's
           no point in warning the user in that case. In any other case
           the user might be able to do something about it. */
#ifndef __APPLE__
        md_print_warn(NULL, fplog,
                      "Can not set thread affinities on the current platform. On NUMA systems this\n"
                      "can cause performance degradation. If you think your platform should support\n"
                      "setting affinities, contact the GROMACS developers.");
#endif  /* __APPLE__ */
        return;
    }

    /* threads on this MPI process or TMPI thread */
    if (cr->duty & DUTY_PP)
    {
        nthread_local = gmx_omp_nthreads_get(emntNonbonded);
    }
    else
    {
        nthread_local = gmx_omp_nthreads_get(emntPME);
    }

    /* map the current process to cores */
    thread0_id_node = 0;
    nthread_node    = nthread_local;
#ifdef GMX_MPI
    if (PAR(cr) || MULTISIM(cr))
    {
        /* We need to determine a scan of the thread counts in this
         * compute node.
         */
        MPI_Comm comm_intra;

        MPI_Comm_split(MPI_COMM_WORLD,
                       gmx_physicalnode_id_hash(), cr->rank_intranode,
                       &comm_intra);
        MPI_Scan(&nthread_local, &thread0_id_node, 1, MPI_INT, MPI_SUM, comm_intra);
        /* MPI_Scan is inclusive, but here we need exclusive */
        thread0_id_node -= nthread_local;
        /* Get the total number of threads on this physical node */
        MPI_Allreduce(&nthread_local, &nthread_node, 1, MPI_INT, MPI_SUM, comm_intra);
        MPI_Comm_free(&comm_intra);
    }
#endif

    if (hw_opt->thread_affinity == threadaffAUTO &&
        nthread_node != hwinfo->nthreads_hw_avail)
    {
        if (nthread_node > 1 && nthread_node < hwinfo->nthreads_hw_avail)
        {
            md_print_warn(cr, fplog,
                          "NOTE: The number of threads is not equal to the number of (logical) cores\n"
                          "      and the -pin option is set to auto: will not pin thread to cores.\n"
                          "      This can lead to significant performance degradation.\n"
                          "      Consider using -pin on (and -pinoffset in case you run multiple jobs).\n");
        }

        return;
    }

    // TODO: Assumes that each thread calls function only once.
    // It would be good to assert this is true at least.
    static int thread_id = -1;
    thread_id++;

    offset = 0;
    if (hw_opt->core_pinning_offset != 0)
    {
        offset = hw_opt->core_pinning_offset;
        if (thread_id == 0)
        {
            md_print_info(cr, fplog, "Applying core pinning offset %d\n", offset);
        }
    }

    rc = get_thread_affinity_layout(fplog, cr, hwinfo,
                                    nthread_node,
                                    offset, &hw_opt->core_pinning_stride,
                                    &locality_order);

    if (rc != 0)
    {
        /* Incompatible layout, don't pin, warning was already issued */
        return;
    }

    int thread_id_node = thread0_id_node + thread_id;
    int index = offset + thread_id_node*hw_opt->core_pinning_stride;
    int core;

    if (locality_order != NULL)
    {
        core = locality_order[index];
    }
    else
    {
        core = index;
    }

    int setaffinity_return = tMPI_Thread_setaffinity_single(tMPI_Thread_self(), core);
    if (setaffinity_return != 0)
    {
        md_print_warn(cr, fplog, "%s", "Setting affinity for thread %d failed\n", thread_id);
    }

    if (debug)
    {
        fprintf(debug, "On rank %2d, thread %2d, index %2d, core %2d the affinity setting returned %d\n",
                cr->nodeid, thread_id, index, core, setaffinity_return);
    }
}

// TBB callback infrastructure for setting individual thread affinities
// Calls gmx_set_thread_affinity to actually set affinity value
class affinity_observer: public tbb::task_scheduler_observer {
    FILE *fplog;
    const t_commrec *cr;
    gmx_hw_opt_t *hw_opt;
    const gmx_hw_info_t *hwinfo;

public:
    affinity_observer(FILE                *f,
                      const t_commrec     *c,
                      gmx_hw_opt_t        *opt,
                      const gmx_hw_info_t *info)
    : fplog(f), cr(c), hw_opt(opt), hwinfo(info) {}

    void on_scheduler_entry( bool )
    {
         gmx_set_thread_affinity(fplog, cr, hw_opt, hwinfo);
    }
};

static affinity_observer *tbb_affinity_observer;

void
gmx_init_tbb_threads(tbb::task_scheduler_init &scheduler,
                          FILE *fplog,
                          const t_commrec     *cr,
                          gmx_hw_opt_t        *hw_opt,
                          const gmx_hw_info_t *hwinfo)
{
    scheduler.initialize(gmx_omp_nthreads_get(emntNonbonded));
    tbb_affinity_observer = new affinity_observer(fplog, cr, hw_opt, hwinfo);
    tbb_affinity_observer->observe(1);
}

/* Check the process affinity mask and if it is found to be non-zero,
 * will honor it and disable mdrun internal affinity setting.
 * Note that this will only work on Linux as we use a GNU feature.
 */
void
gmx_check_thread_affinity_set(FILE            *fplog,
                              const t_commrec *cr,
                              gmx_hw_opt_t    *hw_opt,
                              int  gmx_unused  nthreads_hw_avail,
                              gmx_bool         bAfterOpenmpInit)
{
#ifdef HAVE_SCHED_AFFINITY
    cpu_set_t mask_current;
    int       i, ret, cpu_count, cpu_set;
    gmx_bool  bAllSet;
#endif
#ifdef GMX_LIB_MPI
    gmx_bool  bAllSet_All;
#endif

    assert(hw_opt);
    if (!bAfterOpenmpInit)
    {
        /* Check for externally set OpenMP affinity and turn off internal
         * pinning if any is found. We need to do this check early to tell
         * thread-MPI whether it should do pinning when spawning threads.
         * TODO: the above no longer holds, we should move these checks later
         */
        if (hw_opt->thread_affinity != threadaffOFF)
        {
            char *message;
            if (!gmx_omp_check_thread_affinity(&message))
            {
                /* TODO: with -pin auto we should only warn when using all cores */
                md_print_warn(cr, fplog, "%s", message);
                sfree(message);
                hw_opt->thread_affinity = threadaffOFF;
            }
        }

        /* With thread-MPI this is needed as pinning might get turned off,
         * which needs to be known before starting thread-MPI.
         * With thread-MPI hw_opt is processed here on the master rank
         * and passed to the other ranks later, so we only do this on master.
         */
        if (!SIMMASTER(cr))
        {
            return;
        }
#ifndef GMX_THREAD_MPI
        return;
#endif
    }

#ifdef HAVE_SCHED_GETAFFINITY
    if (hw_opt->thread_affinity == threadaffOFF)
    {
        /* internal affinity setting is off, don't bother checking process affinity */
        return;
    }

    CPU_ZERO(&mask_current);
    if ((ret = sched_getaffinity(0, sizeof(cpu_set_t), &mask_current)) != 0)
    {
        /* failed to query affinity mask, will just return */
        if (debug)
        {
            fprintf(debug, "Failed to query affinity mask (error %d)", ret);
        }
        return;
    }

    /* Before proceeding with the actual check, make sure that the number of
     * detected CPUs is >= the CPUs in the current set.
     * We need to check for CPU_COUNT as it was added only in glibc 2.6. */
#ifdef CPU_COUNT
    if (nthreads_hw_avail < CPU_COUNT(&mask_current))
    {
        if (debug)
        {
            fprintf(debug, "%d hardware threads detected, but %d was returned by CPU_COUNT",
                    nthreads_hw_avail, CPU_COUNT(&mask_current));
        }
        return;
    }
#endif /* CPU_COUNT */

    bAllSet = TRUE;
    for (i = 0; (i < nthreads_hw_avail && i < CPU_SETSIZE); i++)
    {
        bAllSet = bAllSet && (CPU_ISSET(i, &mask_current) != 0);
    }

#ifdef GMX_LIB_MPI
    MPI_Allreduce(&bAllSet, &bAllSet_All, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    bAllSet = bAllSet_All;
#endif

    if (!bAllSet)
    {
        if (hw_opt->thread_affinity == threadaffAUTO)
        {
            if (!bAfterOpenmpInit)
            {
                md_print_warn(cr, fplog,
                              "Non-default thread affinity set, disabling internal thread affinity");
            }
            else
            {
                md_print_warn(cr, fplog,
                              "Non-default thread affinity set probably by the OpenMP library,\n"
                              "disabling internal thread affinity");
            }
            hw_opt->thread_affinity = threadaffOFF;
        }
        else
        {
            /* Only warn once, at the last check (bAfterOpenmpInit==TRUE) */
            if (bAfterOpenmpInit)
            {
                md_print_warn(cr, fplog,
                              "Overriding thread affinity set outside %s\n",
                              ShortProgram());
            }
        }

        if (debug)
        {
            fprintf(debug, "Non-default affinity mask found\n");
        }
    }
    else
    {
        if (debug)
        {
            fprintf(debug, "Default affinity mask found\n");
        }
    }
#endif /* HAVE_SCHED_AFFINITY */
}
