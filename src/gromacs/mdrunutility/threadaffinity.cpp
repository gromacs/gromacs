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

#include "threadaffinity.h"

#include "config.h"

#include <cerrno>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#if HAVE_SCHED_AFFINITY
#    include <sched.h>
#endif

#include "thread_mpi/threads.h"

#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/unique_cptr.h"

namespace
{

class DefaultThreadAffinityAccess : public gmx::IThreadAffinityAccess
{
public:
    bool isThreadAffinitySupported() const override
    {
        return tMPI_Thread_setaffinity_support() == TMPI_SETAFFINITY_SUPPORT_YES;
    }
    bool setCurrentThreadAffinityToCore(int core) override
    {
        const int ret = tMPI_Thread_setaffinity_single(tMPI_Thread_self(), core);
        return ret == 0;
    }
};

//! Global instance of DefaultThreadAffinityAccess
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DefaultThreadAffinityAccess g_defaultAffinityAccess;

} // namespace

gmx::IThreadAffinityAccess::~IThreadAffinityAccess() {}

static bool invalidWithinSimulation(const t_commrec* cr, bool invalidLocally)
{
#if GMX_MPI
    if (cr->nnodes > 1)
    {
        int value = invalidLocally ? 1 : 0;
        int globalValue;
        MPI_Reduce(&value, &globalValue, 1, MPI_INT, MPI_LOR, MAINRANK(cr), cr->mpi_comm_mysim);
        return SIMMAIN(cr) ? (globalValue != 0) : invalidLocally;
    }
#else
    GMX_UNUSED_VALUE(cr);
#endif
    return invalidLocally;
}

static bool get_thread_affinity_layout(const gmx::MDLogger&         mdlog,
                                       const t_commrec*             cr,
                                       const gmx::HardwareTopology& hwTop,
                                       int                          threads,
                                       bool  affinityIsAutoAndNumThreadsIsNotAuto,
                                       int   pin_offset,
                                       int*  pin_stride,
                                       int** localityOrder,
                                       bool* issuedWarning)
{
    bool        bPickPinStride;
    bool        haveTopology;
    bool        invalidValue;
    int         hwMaxThreads  = hwTop.maxThreads();
    std::size_t maxSmtPerCore = 1;

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
        snew(*localityOrder, hwTop.machine().logicalProcessors.size());
        int i = 0;
        for (const auto& s : hwTop.machine().packages)
        {
            for (const auto& c : s.cores)
            {
                maxSmtPerCore = std::max(maxSmtPerCore, c.processingUnits.size());
                for (const auto& pu : c.processingUnits)
                {
                    (*localityOrder)[i++] = pu.osId;
                }
            }
        }
    }
    else
    {
        /* topology information not available or invalid, ignore it */
        *localityOrder = nullptr;
    }
    // Only warn about the first problem per node.  Otherwise, the first test
    // failing would essentially always cause also the other problems get
    // reported, leading to bogus warnings.  The order in the conditionals
    // with this variable is important, since the MPI_Reduce() in
    // invalidWithinSimulation() needs to always happen.
    bool alreadyWarned = false;
    invalidValue       = (hwMaxThreads <= 0);
    if (invalidWithinSimulation(cr, invalidValue))
    {
        /* We don't know anything about the hardware, don't pin */
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "NOTE: No information on available logical cpus, thread pinning disabled.");
        alreadyWarned = true;
    }
    bool validLayout = !invalidValue;

    if (haveTopology)
    {
        invalidValue = (hwMaxThreads < static_cast<int>(hwTop.machine().logicalProcessors.size()));
        if (invalidWithinSimulation(cr, invalidValue) && !alreadyWarned)
        {
            // Don't pin on things that look like containers with shared resources
            // where we are limited to only using a fraction of them
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText(
                            "NOTE: OS CPU limit is lower than logical cpu count, thread pinning "
                            "disabled.");
            alreadyWarned = true;
        }
        validLayout = validLayout && !invalidValue;
    }

    if (affinityIsAutoAndNumThreadsIsNotAuto)
    {
        invalidValue = (threads != hwMaxThreads);
        bool warn    = (invalidValue && threads > 1 && threads < hwMaxThreads);
        if (invalidWithinSimulation(cr, warn) && !alreadyWarned)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText(
                            "NOTE: The number of threads is not equal to the number of (logical) "
                            "cpus\n"
                            "      and the -pin option is set to auto: will not pin threads to "
                            "cpus.\n"
                            "      This can lead to significant performance degradation.\n"
                            "      Consider using -pin on (and -pinoffset in case you run multiple "
                            "jobs).");
            alreadyWarned = true;
        }
        validLayout = validLayout && !invalidValue;
    }

    invalidValue = (threads > hwMaxThreads);
    if (invalidWithinSimulation(cr, invalidValue) && !alreadyWarned)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText("NOTE: Oversubscribing available/permitted CPUs, will not pin threads");
        alreadyWarned = true;
    }
    validLayout = validLayout && !invalidValue;

    invalidValue = (pin_offset + threads > hwMaxThreads);
    if (invalidWithinSimulation(cr, invalidValue) && !alreadyWarned)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "WARNING: Requested offset too large for available logical cpus, thread "
                        "pinning "
                        "disabled.");
        alreadyWarned = true;
    }
    validLayout = validLayout && !invalidValue;

    invalidValue = false;
    /* do we need to choose the pinning stride? */
    bPickPinStride = (*pin_stride == 0);

    if (bPickPinStride)
    {
        // Strides get REALLY yucky on modern hybrid CPUs that might combine
        // performance cores having SMT with efficiency ones that don't.
        // For now we simply start by testing if we can stride it with the maxSmt
        if (haveTopology && pin_offset + threads * static_cast<int>(maxSmtPerCore) <= hwMaxThreads)
        {
            // If all cores are uniform, this will get place one thread per core.
            // If they are not, we hope the performance cores come first, which
            // should get us one thread per those cores at least - and then we
            // might waste some efficiency cores.
            *pin_stride = maxSmtPerCore;
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
            *pin_stride = (hwMaxThreads - pin_offset) / threads;
        }
    }
    else
    {
        /* Check the placement of the thread with the largest index to make sure
         * that the offset & stride doesn't cause pinning beyond the last hardware thread. */
        invalidValue = (pin_offset + (threads - 1) * (*pin_stride) >= hwMaxThreads);
    }
    if (invalidWithinSimulation(cr, invalidValue) && !alreadyWarned)
    {
        /* We are oversubscribing, don't pin */
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "WARNING: Requested stride too large for available logical cpus, thread "
                        "pinning "
                        "disabled.");
        alreadyWarned = true;
    }
    validLayout = validLayout && !invalidValue;

    if (validLayout)
    {
        GMX_LOG(mdlog.info)
                .appendTextFormatted("Pinning threads with a%s logical cpu stride of %d",
                                     bPickPinStride ? "n auto-selected" : " user-specified",
                                     *pin_stride);
    }

    *issuedWarning = alreadyWarned;

    return validLayout;
}

static bool set_affinity(const t_commrec*            cr,
                         int                         nthread_local,
                         int                         intraNodeThreadOffset,
                         int                         offset,
                         int                         core_pinning_stride,
                         const int*                  localityOrder,
                         gmx::IThreadAffinityAccess* affinityAccess)
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
#pragma omp parallel num_threads(nthread_local) reduction(+ : nth_affinity_set)
    {
        try
        {
            int thread_id, thread_id_node;
            int index, core;

            thread_id      = gmx_omp_get_thread_num();
            thread_id_node = intraNodeThreadOffset + thread_id;
            index          = offset + thread_id_node * core_pinning_stride;
            if (localityOrder != nullptr)
            {
                core = localityOrder[index];
            }
            else
            {
                core = index;
            }

            const bool ret = affinityAccess->setCurrentThreadAffinityToCore(core);

            /* store the per-thread success-values of the setaffinity */
            nth_affinity_set += (ret ? 1 : 0);

            if (debug)
            {
                fprintf(debug,
                        "On rank %2d, thread %2d, index %2d, core %2d the affinity setting "
                        "returned %d\n",
                        cr->nodeid,
                        gmx_omp_get_thread_num(),
                        index,
                        core,
                        ret ? 1 : 0);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    if (nth_affinity_set > nthread_local)
    {
        char msg[STRLEN];

        sprintf(msg,
                "Looks like we have set affinity for more threads than "
                "we have (%d > %d)!\n",
                nth_affinity_set,
                nthread_local);
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
#    if GMX_THREAD_MPI
            sprintf(sbuf1, "In tMPI thread #%d: ", cr->nodeid);
#    else /* GMX_LIB_MPI */
            sprintf(sbuf1, "In MPI process #%d: ", cr->nodeid);
#    endif
#endif /* GMX_MPI */
        }

        if (nthread_local > 1)
        {
            sprintf(sbuf2,
                    "for %d/%d thread%s ",
                    nthread_local - nth_affinity_set,
                    nthread_local,
                    nthread_local > 1 ? "s" : "");
        }

        // TODO: This output should also go through mdlog.
        fprintf(stderr, "NOTE: %sAffinity setting %sfailed.\n", sbuf1, sbuf2);
    }
    return allAffinitiesSet;
}

void analyzeThreadsOnThisNode(const gmx::PhysicalNodeCommunicator& physicalNodeComm,
                              int                                  numThreadsOnThisRank,
                              int*                                 numThreadsOnThisNode,
                              int*                                 intraNodeThreadOffset)
{
    *intraNodeThreadOffset = 0;
    *numThreadsOnThisNode  = numThreadsOnThisRank;
#if GMX_MPI
    if (physicalNodeComm.size_ > 1)
    {
        /* We need to determine a scan of the thread counts in this
         * compute node. */
        MPI_Scan(&numThreadsOnThisRank, intraNodeThreadOffset, 1, MPI_INT, MPI_SUM, physicalNodeComm.comm_);
        /* MPI_Scan is inclusive, but here we need exclusive */
        *intraNodeThreadOffset -= numThreadsOnThisRank;
        /* Get the total number of threads on this physical node */
        MPI_Allreduce(
                &numThreadsOnThisRank, numThreadsOnThisNode, 1, MPI_INT, MPI_SUM, physicalNodeComm.comm_);
    }
#else
    GMX_UNUSED_VALUE(physicalNodeComm);
#endif
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
void gmx_set_thread_affinity(const gmx::MDLogger&         mdlog,
                             const t_commrec*             cr,
                             const gmx_hw_opt_t*          hw_opt,
                             const gmx::HardwareTopology& hwTop,
                             int                          numThreadsOnThisRank,
                             int                          numThreadsOnThisNode,
                             int                          intraNodeThreadOffset,
                             gmx::IThreadAffinityAccess*  affinityAccess)
{
    int* localityOrder = nullptr;

    if (hw_opt->threadAffinity == ThreadAffinity::Off)
    {
        /* Nothing to do */
        return;
    }

    if (affinityAccess == nullptr)
    {
        affinityAccess = &g_defaultAffinityAccess;
    }

    /* If the tMPI thread affinity setting is not supported encourage the user
     * to report it as it's either a bug or an exotic platform which we might
     * want to support. */
    if (!affinityAccess->isThreadAffinitySupported())
    {
        /* we know Mac OS does not support setting thread affinity, so there's
           no point in warning the user in that case. In any other case
           the user might be able to do something about it. */
#if !defined(__APPLE__)
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText("NOTE: Cannot set thread affinities on the current platform.");
#endif /* __APPLE__ */
        return;
    }

    int offset              = hw_opt->core_pinning_offset;
    int core_pinning_stride = hw_opt->core_pinning_stride;
    if (offset != 0)
    {
        GMX_LOG(mdlog.warning).appendTextFormatted("Applying core pinning offset %d", offset);
    }

    bool affinityIsAutoAndNumThreadsIsNotAuto =
            (hw_opt->threadAffinity == ThreadAffinity::Auto && !hw_opt->totNumThreadsIsAuto);
    bool                   issuedWarning;
    bool                   validLayout = get_thread_affinity_layout(mdlog,
                                                  cr,
                                                  hwTop,
                                                  numThreadsOnThisNode,
                                                  affinityIsAutoAndNumThreadsIsNotAuto,
                                                  offset,
                                                  &core_pinning_stride,
                                                  &localityOrder,
                                                  &issuedWarning);
    const gmx::sfree_guard localityOrderGuard(localityOrder);

    bool allAffinitiesSet;
    if (validLayout)
    {
        allAffinitiesSet = set_affinity(
                cr, numThreadsOnThisRank, intraNodeThreadOffset, offset, core_pinning_stride, localityOrder, affinityAccess);
    }
    else
    {
        // Produce the warning if any rank fails.
        allAffinitiesSet = false;
    }
    if (invalidWithinSimulation(cr, !allAffinitiesSet) && !issuedWarning)
    {
        GMX_LOG(mdlog.warning).asParagraph().appendText("NOTE: Thread affinity was not set.");
    }
}

/* Detects and returns whether we have the default affinity mask
 *
 * Returns true when we can query thread affinities and CPU count is
 * consistent and we have default affinity mask on all ranks.
 *
 * Should be called simultaneously by all MPI ranks.
 */
static bool detectDefaultAffinityMask(const int maxThreads, MPI_Comm world)
{
    bool detectedDefaultAffinityMask = true;

#if HAVE_SCHED_AFFINITY
    cpu_set_t mask_current;
    CPU_ZERO(&mask_current);
    int ret;
    if ((ret = sched_getaffinity(0, sizeof(cpu_set_t), &mask_current)) != 0)
    {
        /* failed to query affinity mask, will just return */
        if (debug)
        {
            fprintf(debug, "Failed to query affinity mask (error %d)", ret);
        }
        detectedDefaultAffinityMask = false;
    }

    if (detectedDefaultAffinityMask)
    {
        /* Here we check whether all bits of the affinity mask are set.
         * Note that this mask can change while the program is executing,
         * e.g. by the system dynamically enabling or disabling cores or
         * by some scheduler changing the affinities. Thus we can not assume
         * that the result is the same on all ranks within a node
         * when running this code at different times.
         */
        bool allBitsAreSet = true;
        for (int i = 0; (i < maxThreads && i < CPU_SETSIZE); i++)
        {
            allBitsAreSet = allBitsAreSet && (CPU_ISSET(i, &mask_current) != 0);
        }
        if (debug)
        {
            fprintf(debug, "%s affinity mask found\n", allBitsAreSet ? "Default" : "Non-default");
        }
        if (!allBitsAreSet)
        {
            detectedDefaultAffinityMask = false;
        }
    }
#else
    GMX_UNUSED_VALUE(maxThreads);
#endif /* HAVE_SCHED_AFFINITY */

#if GMX_MPI
    int mpiIsInitialized;
    MPI_Initialized(&mpiIsInitialized);
    if (mpiIsInitialized)
    {
        bool maskToReduce = detectedDefaultAffinityMask;
        MPI_Allreduce(&maskToReduce, &detectedDefaultAffinityMask, 1, MPI_C_BOOL, MPI_LAND, world);
    }
#else
    GMX_UNUSED_VALUE(world);
#endif

    return detectedDefaultAffinityMask;
}

/* Check the process affinity mask and if it is found to be non-zero,
 * will honor it and disable mdrun internal affinity setting.
 * Note that this will only work on Linux as we use a GNU feature.
 */
void gmx_check_thread_affinity_set(const gmx::MDLogger& mdlog,
                                   gmx_hw_opt_t*        hw_opt,
                                   int                  ncpus,
                                   gmx_bool             bAfterOpenmpInit,
                                   MPI_Comm             world)
{
    GMX_RELEASE_ASSERT(hw_opt, "hw_opt must be a non-NULL pointer");

    if (!bAfterOpenmpInit)
    {
        /* Check for externally set OpenMP affinity and turn off internal
         * pinning if any is found. We need to do this check early to tell
         * thread-MPI whether it should do pinning when spawning threads.
         * TODO: the above no longer holds, we should move these checks later
         */
        if (hw_opt->threadAffinity != ThreadAffinity::Off)
        {
            char* message;
            if (!gmx_omp_check_thread_affinity(&message))
            {
                /* We only pin automatically with totNumThreadsIsAuto=true */
                if (hw_opt->threadAffinity == ThreadAffinity::On || hw_opt->totNumThreadsIsAuto)
                {
                    GMX_LOG(mdlog.warning).asParagraph().appendText(message);
                }
                sfree(message);
                hw_opt->threadAffinity = ThreadAffinity::Off;
            }
        }
    }

    if (!detectDefaultAffinityMask(ncpus, world))
    {
        if (hw_opt->threadAffinity == ThreadAffinity::Auto)
        {
            if (!bAfterOpenmpInit)
            {
                GMX_LOG(mdlog.warning)
                        .asParagraph()
                        .appendText(
                                "Non-default thread affinity set, disabling internal thread "
                                "affinity");
            }
            else
            {
                GMX_LOG(mdlog.warning)
                        .asParagraph()
                        .appendText(
                                "Non-default thread affinity set probably by the OpenMP library,\n"
                                "disabling internal thread affinity");
            }
            hw_opt->threadAffinity = ThreadAffinity::Off;
        }
        else
        {
            /* Only warn once, at the last check (bAfterOpenmpInit==TRUE) */
            if (bAfterOpenmpInit)
            {
                GMX_LOG(mdlog.warning)
                        .asParagraph()
                        .appendTextFormatted("Overriding thread affinity set outside %s",
                                             gmx::getProgramContext().displayName());
            }
        }
    }
}
