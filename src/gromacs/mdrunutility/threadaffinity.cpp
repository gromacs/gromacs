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

#include "threadaffinity.h"

#include "config.h"

#include <cerrno>
#include <cstdio>
#include <cstring>

#if HAVE_SCHED_AFFINITY
#  include <sched.h>
#  include <sys/syscall.h>
#endif

#include "thread_mpi/threads.h"

#include "gromacs/gmxlib/md_logging.h"
#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/scoped_cptr.h"
#include "gromacs/utility/smalloc.h"


static bool invalidWithinSimulation(const t_commrec *cr, bool invalidLocally)
{
#if GMX_MPI
    if (cr->nnodes > 1)
    {
        int value = invalidLocally ? 1 : 0;
        int globalValue;
        MPI_Reduce(&value, &globalValue, 1, MPI_INT, MPI_LOR, MASTERRANK(cr),
                   cr->mpi_comm_mysim);
        return SIMMASTER(cr) ? (globalValue != 0) : invalidLocally;
    }
#else
    GMX_UNUSED_VALUE(cr);
#endif
    return invalidLocally;
}

static bool
get_thread_affinity_layout(FILE *fplog,
                           const t_commrec *cr,
                           const gmx_hw_info_t * hwinfo,
                           int   threads,
                           int pin_offset, int * pin_stride,
                           int **localityOrder)
{
    int                          hwThreads;
    int                          hwThreadsPerCore = 0;
    bool                         bPickPinStride;
    bool                         haveTopology;
    bool                         invalidValue;

    const gmx::HardwareTopology &hwTop = *hwinfo->hardwareTopology;

    haveTopology = (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Basic);

    if (pin_offset < 0)
    {
        gmx_fatal(FARGS, "Negative thread pinning offset requested");
    }
    if (*pin_stride < 0)
    {
        gmx_fatal(FARGS, "Negative thread pinning stride requested");
    }

    if (haveTopology)
    {
        hwThreads           = hwTop.machine().logicalProcessorCount;
        // Just use the value for the first core
        hwThreadsPerCore    = hwTop.machine().sockets[0].cores[0].hwThreads.size();
        snew(*localityOrder, hwThreads);
        int i = 0;
        for (auto &s : hwTop.machine().sockets)
        {
            for (auto &c : s.cores)
            {
                for (auto &t : c.hwThreads)
                {
                    (*localityOrder)[i++] = t.logicalProcessorId;
                }
            }
        }
    }
    else
    {
        /* topology information not available or invalid, ignore it */
        hwThreads       = hwinfo->nthreads_hw_avail;
        *localityOrder  = NULL;
    }
    bool validLayout = (hwThreads > 0);
    if (!validLayout)
    {
        /* We don't know anything about the hardware, don't pin */
        md_print_warn(cr, fplog,
                      "NOTE: No information on available cores, thread pinning disabled.\n");
    }

    invalidValue = (threads > hwThreads);
    // Only warn about the first problem per node.  Otherwise, the first test
    // failing would essentially always cause also the other problems get
    // reported, leading to bogus warnings.  The order in the conditional is
    // important, since the MPI_Reduce() needs to always happen.
    if (invalidWithinSimulation(cr, invalidValue) && validLayout)
    {
        md_print_warn(cr, fplog,
                      "NOTE: Oversubscribing a CPU, will not pin threads.\n");
    }
    validLayout = validLayout && !invalidValue;

    invalidValue = (pin_offset + threads > hwThreads);
    if (invalidWithinSimulation(cr, invalidValue) && validLayout)
    {
        md_print_warn(cr, fplog,
                      "WARNING: Requested offset too large for available cores, thread pinning disabled.\n");

    }
    validLayout = validLayout && !invalidValue;

    invalidValue   = false;
    /* do we need to choose the pinning stride? */
    bPickPinStride = (*pin_stride == 0);

    if (bPickPinStride)
    {
        if (haveTopology && pin_offset + threads*hwThreadsPerCore <= hwThreads)
        {
            /* Put one thread on each physical core */
            *pin_stride = hwThreadsPerCore;
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
            *pin_stride = (hwThreads - pin_offset)/threads;
        }
    }
    else
    {
        /* Check the placement of the thread with the largest index to make sure
         * that the offset & stride doesn't cause pinning beyond the last hardware thread. */
        invalidValue = (pin_offset + (threads-1)*(*pin_stride) >= hwThreads);
    }
    if (invalidWithinSimulation(cr, invalidValue) && validLayout)
    {
        /* We are oversubscribing, don't pin */
        md_print_warn(cr, fplog,
                      "WARNING: Requested stride too large for available cores, thread pinning disabled.\n");

    }
    validLayout = validLayout && !invalidValue;

    if (validLayout && fplog != NULL)
    {
        fprintf(fplog, "Pinning threads with a%s logical core stride of %d\n",
                bPickPinStride ? "n auto-selected" : " user-specified",
                *pin_stride);
    }

    return validLayout;
}

static bool set_affinity(const t_commrec *cr, int nthread_local, int thread0_id_node,
                         int offset, int core_pinning_stride, int *localityOrder)
{
    // Set the per-thread affinity. In order to be able to check the success
    // of affinity settings, we will set nth_affinity_set to 1 on threads
    // where the affinity setting succeded and to 0 where it failed.
    // Reducing these 0/1 values over the threads will give the total number
    // of threads on which we succeeded.

    // To avoid warnings from the static analyzer we initialize nth_affinity_set
    // to zero outside the OpenMP block, and then add to it inside the block.
    // The value will still always be 0 or 1 from each thread.
    int nth_affinity_set = 0;
#pragma omp parallel num_threads(nthread_local) reduction(+:nth_affinity_set)
    {
        try
        {
            int      thread_id, thread_id_node;
            int      index, core;
            gmx_bool setaffinity_ret;

            thread_id      = gmx_omp_get_thread_num();
            thread_id_node = thread0_id_node + thread_id;
            index          = offset + thread_id_node*core_pinning_stride;
            if (localityOrder != nullptr)
            {
                core = localityOrder[index];
            }
            else
            {
                core = index;
            }

            setaffinity_ret = tMPI_Thread_setaffinity_single(tMPI_Thread_self(), core);

            /* store the per-thread success-values of the setaffinity */
            nth_affinity_set += (setaffinity_ret == 0);

            if (debug)
            {
                fprintf(debug, "On rank %2d, thread %2d, index %2d, core %2d the affinity setting returned %d\n",
                        cr->nodeid, gmx_omp_get_thread_num(), index, core, setaffinity_ret);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    if (nth_affinity_set > nthread_local)
    {
        char msg[STRLEN];

        sprintf(msg, "Looks like we have set affinity for more threads than "
                "we have (%d > %d)!\n", nth_affinity_set, nthread_local);
        gmx_incons(msg);
    }

    /* check & warn if some threads failed to set their affinities */
    const bool allAffinitiesSet = (nth_affinity_set == nthread_local);
    if (!allAffinitiesSet)
    {
        char sbuf1[STRLEN], sbuf2[STRLEN];

        /* sbuf1 contains rank info, while sbuf2 OpenMP thread info */
        sbuf1[0] = sbuf2[0] = '\0';
        /* Only add rank info if we have more than one rank. */
        if (cr->nnodes > 1)
        {
#if GMX_MPI
#if GMX_THREAD_MPI
            sprintf(sbuf1, "In tMPI thread #%d: ", cr->nodeid);
#else           /* GMX_LIB_MPI */
            sprintf(sbuf1, "In MPI process #%d: ", cr->nodeid);
#endif
#endif          /* GMX_MPI */
        }

        if (nthread_local > 1)
        {
            sprintf(sbuf2, "for %d/%d thread%s ",
                    nthread_local - nth_affinity_set, nthread_local,
                    nthread_local > 1 ? "s" : "");
        }

        fprintf(stderr, "NOTE: %sAffinity setting %sfailed.\n", sbuf1, sbuf2);
    }
    return allAffinitiesSet;
}

/* Set CPU affinity. Can be important for performance.
   On some systems (e.g. Cray) CPU Affinity is set by default.
   But default assigning doesn't work (well) with only some ranks
   having threads. This causes very low performance.
   External tools have cumbersome syntax for setting affinity
   in the case that only some ranks have threads.
   Thus it is important that GROMACS sets the affinity internally
   if only PME is using threads.
 */
void
gmx_set_thread_affinity(FILE                *fplog,
                        const t_commrec     *cr,
                        const gmx_hw_opt_t  *hw_opt,
                        const gmx_hw_info_t *hwinfo)
{
    int        thread0_id_node,
               nthread_local, nthread_node;
    int        offset;
    int *      localityOrder = nullptr;

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
        /* we know Mac OS & BlueGene do not support setting thread affinity, so there's
           no point in warning the user in that case. In any other case
           the user might be able to do something about it. */
#if !defined(__APPLE__) && !defined(__bg__)
        md_print_warn(cr, fplog,
                      "NOTE: Cannot set thread affinities on the current platform.\n");
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
#if GMX_MPI
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

    offset = 0;
    if (hw_opt->core_pinning_offset != 0)
    {
        offset = hw_opt->core_pinning_offset;
        md_print_info(cr, fplog, "Applying core pinning offset %d\n", offset);
    }

    int  core_pinning_stride = hw_opt->core_pinning_stride;
    bool validLayout
        = get_thread_affinity_layout(fplog, cr, hwinfo, nthread_node, offset,
                                     &core_pinning_stride, &localityOrder);
    gmx::scoped_guard_sfree localityOrderGuard(localityOrder);

    bool                    allAffinitiesSet;
    if (validLayout)
    {
        allAffinitiesSet = set_affinity(cr, nthread_local, thread0_id_node,
                                        offset, core_pinning_stride, localityOrder);
    }
    else
    {
        // Produce the warning if any rank fails.
        allAffinitiesSet = false;
    }
    if (invalidWithinSimulation(cr, !allAffinitiesSet))
    {
        md_print_warn(cr, fplog,
                      "NOTE: Thread affinity setting failed. This can cause performance degradation.\n"
                      "      If you think your settings are correct, ask on the gmx-users list.\n");
    }
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
    GMX_RELEASE_ASSERT(hw_opt, "hw_opt must be a non-NULL pointer");

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
#if !GMX_THREAD_MPI
        return;
#endif
    }

#if HAVE_SCHED_AFFINITY
    int       ret;
    cpu_set_t mask_current;

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

    gmx_bool bAllSet = TRUE;
    for (int i = 0; (i < nthreads_hw_avail && i < CPU_SETSIZE); i++)
    {
        bAllSet = bAllSet && (CPU_ISSET(i, &mask_current) != 0);
    }

#if GMX_LIB_MPI
    gmx_bool  bAllSet_All;

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
                              gmx::getProgramContext().displayName());
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
