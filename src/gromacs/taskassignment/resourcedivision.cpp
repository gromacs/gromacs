/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Defines utility functionality for dividing resources and
 * checking for consistency and usefulness.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */

#include "gmxpre.h"

#include "resourcedivision.h"

#include "config.h"

#include <stdlib.h>
#include <string.h>

#include <algorithm>

#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/taskassignment/hardwareassign.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"


/* DISCLAIMER: All the atom count and thread numbers below are heuristic.
 * The real switching points will depend on the system simulation,
 * the algorithms used and the hardware it's running on, as well as if there
 * are other jobs running on the same machine. We try to take into account
 * factors that have a large influence, such as recent Intel CPUs being
 * much better at wide multi-threading. The remaining factors should
 * (hopefully) have a small influence, such that the performance just before
 * and after a switch point doesn't change too much.
 */

//! Constant used to help minimize preprocessed code
static const bool bHasOmpSupport = GMX_OPENMP;

#if GMX_THREAD_MPI
/* The minimum number of atoms per tMPI thread. With fewer atoms than this,
 * the number of threads will get lowered.
 */
static const int min_atoms_per_mpi_thread =  90;
static const int min_atoms_per_gpu        = 900;
#endif /* GMX_THREAD_MPI */

/**@{*/
/*! \brief Constants for implementing default divisions of threads */

/* TODO choose nthreads_omp based on hardware topology
   when we have a hardware topology detection library */
/* First we consider the case of no MPI (1 MPI rank).
 * In general, when running up to 8 threads, OpenMP should be faster.
 * Note: on AMD Bulldozer we should avoid running OpenMP over two dies.
 * On Intel>=Nehalem running OpenMP on a single CPU is always faster,
 * even on two CPUs it's usually faster (but with many OpenMP threads
 * it could be faster not to use HT, currently we always use HT).
 * On Nehalem/Westmere we want to avoid running 16 threads over
 * two CPUs with HT, so we need a limit<16; thus we use 12.
 * A reasonable limit for Intel Sandy and Ivy bridge,
 * not knowing the topology, is 16 threads.
 * Below we check for Intel and AVX, which for now includes
 * Sandy/Ivy Bridge, Has/Broadwell. By checking for AVX instead of
 * model numbers we ensure also future Intel CPUs are covered.
 */
const int nthreads_omp_faster_default   =  8;
const int nthreads_omp_faster_Nehalem   = 12;
const int nthreads_omp_faster_Intel_AVX = 16;
const int nthreads_omp_faster_AMD_Ryzen = 16;
/* For CPU only runs the fastest options are usually MPI or OpenMP only.
 * With one GPU, using MPI only is almost never optimal, so we need to
 * compare running pure OpenMP with combined MPI+OpenMP. This means higher
 * OpenMP threads counts can still be ok. Multiplying the numbers above
 * by a factor of 2 seems to be a good estimate.
 */
const int nthreads_omp_faster_gpu_fac   =  2;

/* This is the case with MPI (2 or more MPI PP ranks).
 * By default we will terminate with a fatal error when more than 8
 * OpenMP thread are (indirectly) requested, since using less threads
 * nearly always results in better performance.
 * With thread-mpi and multiple GPUs or one GPU and too many threads
 * we first try 6 OpenMP threads and then less until the number of MPI ranks
 * is divisible by the number of GPUs.
 */
#if GMX_OPENMP && GMX_MPI
const int nthreads_omp_mpi_ok_max              =  8;
const int nthreads_omp_mpi_ok_min_cpu          =  1;
#endif
const int nthreads_omp_mpi_ok_min_gpu          =  2;
const int nthreads_omp_mpi_target_max          =  6;

/**@}*/

/*! \brief Returns the maximum OpenMP thread count for which using a single MPI rank
 * should be faster than using multiple ranks with the same total thread count.
 */
static int nthreads_omp_faster(const gmx::CpuInfo &cpuInfo, gmx_bool bUseGPU)
{
    int nth;

    if (cpuInfo.vendor() == gmx::CpuInfo::Vendor::Intel &&
        cpuInfo.feature(gmx::CpuInfo::Feature::X86_Avx))
    {
        nth = nthreads_omp_faster_Intel_AVX;
    }
    else if (gmx::cpuIsX86Nehalem(cpuInfo))
    {
        // Intel Nehalem
        nth = nthreads_omp_faster_Nehalem;
    }
    else if (cpuInfo.vendor() == gmx::CpuInfo::Vendor::Amd && cpuInfo.family() >= 23)
    {
        // AMD Ryzen
        nth = nthreads_omp_faster_AMD_Ryzen;
    }
    else
    {
        nth = nthreads_omp_faster_default;
    }

    if (bUseGPU)
    {
        nth *= nthreads_omp_faster_gpu_fac;
    }

    nth = std::min(nth, GMX_OPENMP_MAX_THREADS);

    return nth;
}

/*! \brief Returns that maximum OpenMP thread count that passes the efficiency check */
gmx_unused static int nthreads_omp_efficient_max(int gmx_unused       nrank,
                                                 const gmx::CpuInfo  &cpuInfo,
                                                 gmx_bool             bUseGPU)
{
#if GMX_OPENMP && GMX_MPI
    if (nrank > 1)
    {
        return nthreads_omp_mpi_ok_max;
    }
    else
#endif
    {
        return nthreads_omp_faster(cpuInfo, bUseGPU);
    }
}

/*! \brief Return the number of thread-MPI ranks to use.
 * This is chosen such that we can always obey our own efficiency checks.
 */
gmx_unused static int get_tmpi_omp_thread_division(const gmx_hw_info_t *hwinfo,
                                                   const gmx_hw_opt_t  &hw_opt,
                                                   int                  nthreads_tot,
                                                   int                  ngpu)
{
    int                 nrank;
    const gmx::CpuInfo &cpuInfo = *hwinfo->cpuInfo;

    GMX_RELEASE_ASSERT(nthreads_tot > 0, "There must be at least one thread per rank");

    /* There are no separate PME nodes here, as we ensured in
     * check_and_update_hw_opt that nthreads_tmpi>0 with PME nodes
     * and a conditional ensures we would not have ended up here.
     * Note that separate PME nodes might be switched on later.
     */
    if (ngpu > 0)
    {
        nrank = ngpu;

        /* When the user sets nthreads_omp, we can end up oversubscribing CPU cores
         * if we simply start as many ranks as GPUs. To avoid this, we start as few
         * tMPI ranks as necessary to avoid oversubscription and instead leave GPUs idle.
         * If the user does not set the number of OpenMP threads, nthreads_omp==0 and
         * this code has no effect.
         */
        GMX_RELEASE_ASSERT(hw_opt.nthreads_omp >= 0, "nthreads_omp is negative, but previous checks should have prevented this");
        while (nrank*hw_opt.nthreads_omp > hwinfo->nthreads_hw_avail && nrank > 1)
        {
            nrank--;
        }

        if (nthreads_tot < nrank)
        {
            /* #thread < #gpu is very unlikely, but if so: waste gpu(s) */
            nrank = nthreads_tot;
        }
        else if (nthreads_tot > nthreads_omp_faster(cpuInfo, ngpu > 0) ||
                 (ngpu > 1 && nthreads_tot/ngpu > nthreads_omp_mpi_target_max))
        {
            /* The high OpenMP thread count will likely result in sub-optimal
             * performance. Increase the rank count to reduce the thread count
             * per rank. This will lead to GPU sharing by MPI ranks/threads.
             */
            int nshare;

            /* Increase the rank count as long as have we more than 6 OpenMP
             * threads per rank or the number of hardware threads is not
             * divisible by the rank count. Don't go below 2 OpenMP threads.
             */
            nshare = 1;
            do
            {
                nshare++;
                nrank = ngpu*nshare;
            }
            while (nthreads_tot/nrank > nthreads_omp_mpi_target_max ||
                   (nthreads_tot/(ngpu*(nshare + 1)) >= nthreads_omp_mpi_ok_min_gpu && nthreads_tot % nrank != 0));
        }
    }
    else if (hw_opt.nthreads_omp > 0)
    {
        /* Here we could oversubscribe, when we do, we issue a warning later */
        nrank = std::max(1, nthreads_tot/hw_opt.nthreads_omp);
    }
    else
    {
        if (nthreads_tot <= nthreads_omp_faster(cpuInfo, ngpu > 0))
        {
            /* Use pure OpenMP parallelization */
            nrank = 1;
        }
        else
        {
            /* Don't use OpenMP parallelization */
            nrank = nthreads_tot;
        }
    }

    return nrank;
}


#if GMX_THREAD_MPI

static bool
gmxSmtIsEnabled(const gmx::HardwareTopology &hwTop)
{
    return (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Basic && hwTop.machine().sockets[0].cores[0].hwThreads.size() > 1);
}

namespace
{

class SingleRankChecker
{
    public:
        //! Constructor
        SingleRankChecker() : value_(false), reasons_() {}
        /*! \brief Call this function for each possible condition
            under which a single rank is required, along with a string
            describing the constraint when it is applied. */
        void applyConstraint(bool condition, const char *description)
        {
            if (condition)
            {
                value_ = true;
                reasons_.push_back(gmx::formatString("%s only supports a single rank.", description));
            }
        }
        //! After applying any conditions, is a single rank required?
        bool mustUseOneRank() const
        {
            return value_;
        }
        /*! \brief Return a formatted string to use when writing a
            message when a single rank is required, (or empty if no
            constraint exists.) */
        std::string getMessage() const
        {
            return formatAndJoin(reasons_, "\n", gmx::IdentityFormatter());
        }
    private:
        bool                     value_;
        std::vector<std::string> reasons_;
};

} // namespace

/* Get the number of MPI ranks to use for thread-MPI based on how many
 * were requested, which algorithms we're using,
 * and how many particles there are.
 * At the point we have already called check_and_update_hw_opt.
 * Thus all options should be internally consistent and consistent
 * with the hardware, except that ntmpi could be larger than #GPU.
 */
int get_nthreads_mpi(const gmx_hw_info_t    *hwinfo,
                     gmx_hw_opt_t           *hw_opt,
                     int                     numPmeRanks,
                     bool                    nonbondedOnGpu,
                     const t_inputrec       *inputrec,
                     const gmx_mtop_t       *mtop,
                     const gmx::MDLogger    &mdlog,
                     bool                    doMembed)
{
    int                          nthreads_hw, nthreads_tot_max, nrank, ngpu;
    int                          min_atoms_per_mpi_rank;

    const gmx::CpuInfo          &cpuInfo = *hwinfo->cpuInfo;
    const gmx::HardwareTopology &hwTop   = *hwinfo->hardwareTopology;

    /* If the user made a GPU task assignment, that sets the number of thread-MPI ranks. */
    auto userGpuTaskAssignment = gmx::parseGpuTaskAssignment(hw_opt->gpuIdTaskAssignment);
    int  numGpuIdsSupplied     = static_cast<int>(userGpuTaskAssignment.size());

    /* TODO Here we handle the case where the user set GPU IDs, and
       further below we handle the case where the algorithm does not
       support multiple ranks. We need also to handle the case where
       the user set multiple GPU IDs for an algorithm that cannot
       handle multiple ranks. */
    if (hw_opt->nthreads_tmpi < 1 && numGpuIdsSupplied > 0)
    {
        /* If the user chose both mdrun -nt -gpu_id, is that consistent? */
        if (numPmeRanks <= 0)
        {
            if (hw_opt->nthreads_tot > 0 &&
                (hw_opt->nthreads_tot % numGpuIdsSupplied) != 0)
            {
                gmx_fatal(FARGS, "Cannot run %d total threads with %d GPU ranks. Choose the total number of threads to be a multiple of the number of GPU ranks.", hw_opt->nthreads_tot, numGpuIdsSupplied);
            }
            return numGpuIdsSupplied;
        }
        else
        {
            gmx_fatal(FARGS, "The combination of choosing a number of PME ranks, and specific GPU IDs "
                      "is not supported. Use also -ntmpi and/or -ntomp and -ntomp_pme to specify what "
                      "distribution of threads to ranks you require.");
        }
    }

    {
        /* Check if an algorithm does not support parallel simulation.  */
        // TODO This might work better if e.g. implemented algorithms
        // had to define a function that returns such requirements,
        // and a description string.
        SingleRankChecker checker;
        checker.applyConstraint(inputrec->eI == eiLBFGS, "L-BFGS minimization");
        checker.applyConstraint(inputrec->coulombtype == eelEWALD, "Plain Ewald electrostatics");
        checker.applyConstraint(doMembed, "Membrane embedding");
        if (checker.mustUseOneRank())
        {
            std::string message = checker.getMessage();
            if (hw_opt->nthreads_tmpi > 1)
            {
                gmx_fatal(FARGS, "%s However, you asked for more than 1 thread-MPI rank, so mdrun cannot continue. Choose a single rank, or a different algorithm.", message.c_str());
            }
            GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted("%s Choosing to use only a single thread-MPI rank.", message.c_str());

            if (numGpuIdsSupplied > 1)
            {
                gmx_fatal(FARGS, "You supplied %d GPU IDs but only 1 rank can be used "
                          "by this simulation. Supply only one GPU ID.", numGpuIdsSupplied);
            }
            return 1;
        }
    }

    if (hw_opt->nthreads_tmpi > 0)
    {
        if (numPmeRanks <= 0)
        {
            int numPpRanks = hw_opt->nthreads_tmpi;
            if ((numGpuIdsSupplied > 0) &&
                (numGpuIdsSupplied != numPpRanks))
            {
                gmx_fatal(FARGS, "Cannot run %d thread-MPI total ranks with %d "
                          "GPU IDs supplied. The number of particle-particle (PP) ranks and the "
                          "number of GPU IDs must match.", hw_opt->nthreads_tmpi, numGpuIdsSupplied);
            }
        }
        else
        {
            int numPpRanks = hw_opt->nthreads_tmpi - numPmeRanks;
            if ((numGpuIdsSupplied > 0) &&
                (numGpuIdsSupplied != numPpRanks))
            {
                gmx_fatal(FARGS, "Cannot run %d thread-MPI total ranks with %d PME ranks and %d "
                          "GPU IDs supplied. The number of particle-particle ranks and the "
                          "number of GPU IDs must match.", hw_opt->nthreads_tmpi, numPmeRanks, numGpuIdsSupplied);
            }
        }
        /* Trivial, return the user's choice right away */
        return hw_opt->nthreads_tmpi;
    }
    GMX_RELEASE_ASSERT(numGpuIdsSupplied == 0,
                       "If mdrun -gpu_id had information, the number of ranks should have already been chosen");

    // Now implement automatic selection of number of thread-MPI ranks
    nthreads_hw = hwinfo->nthreads_hw_avail;

    if (nthreads_hw <= 0)
    {
        /* This should normally not happen, but if it does, we handle it */
        gmx_fatal(FARGS, "The number of available hardware threads can not be detected, please specify the number of MPI ranks and the number of OpenMP threads (if supported) manually with options -ntmpi and -ntomp, respectively");
    }

    /* How many total (#tMPI*#OpenMP) threads can we start? */
    if (hw_opt->nthreads_tot > 0)
    {
        nthreads_tot_max = hw_opt->nthreads_tot;
    }
    else
    {
        nthreads_tot_max = nthreads_hw;
    }

    /* nonbondedOnGpu might be false e.g. because this simulation uses
     * the group scheme, or is a rerun with energy groups. */
    ngpu = (nonbondedOnGpu ? hwinfo->gpu_info.n_dev_compatible : 0);

    if (inputrec->cutoff_scheme == ecutsGROUP)
    {
        /* We checked this before, but it doesn't hurt to do it once more */
        GMX_RELEASE_ASSERT(hw_opt->nthreads_omp == 1, "The group scheme only supports one OpenMP thread per rank");
    }

    nrank =
        get_tmpi_omp_thread_division(hwinfo, *hw_opt, nthreads_tot_max, ngpu);

    if (inputrec->eI == eiNM || EI_TPI(inputrec->eI))
    {
        /* Dims/steps are divided over the nodes iso splitting the atoms.
         * With NM we can't have more ranks than #atoms*#dim. With TPI it's
         * unlikely we have fewer atoms than ranks, and if so, communication
         * would become a bottleneck, so we set the limit to 1 atom/rank.
         */
        min_atoms_per_mpi_rank = 1;
    }
    else
    {
        if (ngpu >= 1)
        {
            min_atoms_per_mpi_rank = min_atoms_per_gpu;
        }
        else
        {
            min_atoms_per_mpi_rank = min_atoms_per_mpi_thread;
        }
    }

    if (mtop->natoms/nrank < min_atoms_per_mpi_rank)
    {
        int nrank_new;

        /* the rank number was chosen automatically, but there are too few
           atoms per rank, so we need to reduce the rank count */
        nrank_new = std::max(1, mtop->natoms/min_atoms_per_mpi_rank);

        /* Avoid partial use of Hyper-Threading */
        if (gmxSmtIsEnabled(hwTop) &&
            nrank_new > nthreads_hw/2 && nrank_new < nthreads_hw)
        {
            nrank_new = nthreads_hw/2;
        }

        /* If the user specified the total thread count, ensure this is
         * divisible by the number of ranks.
         * It is quite likely that we have too many total threads compared
         * to the size of the system, but if the user asked for this many
         * threads we should respect that.
         */
        while (hw_opt->nthreads_tot > 0 &&
               hw_opt->nthreads_tot % nrank_new != 0)
        {
            nrank_new--;
        }

        /* Avoid large prime numbers in the rank count */
        if (nrank_new >= 6)
        {
            /* Use only 6,8,10 with additional factors of 2 */
            int fac;

            fac = 2;
            while (3*fac*2 <= nrank_new)
            {
                fac *= 2;
            }

            nrank_new = (nrank_new/fac)*fac;
        }
        else
        {
            /* Avoid 5, since small system won't fit 5 domains along
             * a dimension. This might lead to waisting some cores, but this
             * will have a small impact in this regime of very small systems.
             */
            if (nrank_new == 5)
            {
                nrank_new = 4;
            }
        }

        nrank = nrank_new;

        /* We reduced the number of tMPI ranks, which means we might violate
         * our own efficiency checks if we simply use all hardware threads.
         */
        if (bHasOmpSupport && hw_opt->nthreads_omp <= 0 && hw_opt->nthreads_tot <= 0)
        {
            /* The user set neither the total nor the OpenMP thread count,
             * we should use all hardware threads, unless we will violate
             * our own efficiency limitation on the thread count.
             */
            int  nt_omp_max;

            nt_omp_max = nthreads_omp_efficient_max(nrank, cpuInfo, ngpu >= 1);

            if (nrank*nt_omp_max < hwinfo->nthreads_hw_avail)
            {
                /* Limit the number of OpenMP threads to start */
                hw_opt->nthreads_omp = nt_omp_max;
            }
        }

        fprintf(stderr, "\n");
        fprintf(stderr, "NOTE: Parallelization is limited by the small number of atoms,\n");
        fprintf(stderr, "      only starting %d thread-MPI ranks.\n", nrank);
        fprintf(stderr, "      You can use the -nt and/or -ntmpi option to optimize the number of threads.\n\n");
    }

    return nrank;
}
#endif /* GMX_THREAD_MPI */


void check_resource_division_efficiency(const gmx_hw_info_t *hwinfo,
                                        int                  numTotalThreads,
                                        bool                 willUsePhysicalGpu,
                                        gmx_bool             bNtOmpOptionSet,
                                        t_commrec           *cr,
                                        const gmx::MDLogger &mdlog)
{
#if GMX_OPENMP && GMX_MPI
    int         nth_omp_min, nth_omp_max;
    char        buf[1000];
#if GMX_THREAD_MPI
    const char *mpi_option = " (option -ntmpi)";
#else
    const char *mpi_option = "";
#endif

    /* This function should be called after thread-MPI (when configured) and
     * OpenMP have been initialized. Check that here.
     */
#if GMX_THREAD_MPI
    GMX_RELEASE_ASSERT(nthreads_omp_faster_default >= nthreads_omp_mpi_ok_max, "Inconsistent OpenMP thread count default values");
#endif
    GMX_RELEASE_ASSERT(gmx_omp_nthreads_get(emntDefault) >= 1, "Must have at least one OpenMP thread");

    nth_omp_min = gmx_omp_nthreads_get(emntDefault);
    nth_omp_max = gmx_omp_nthreads_get(emntDefault);

    bool anyRankIsUsingGpus = willUsePhysicalGpu;
    /* Thread-MPI seems to have a bug with reduce on 1 node, so use a cond. */
    if (cr->nnodes + cr->npmenodes > 1)
    {
        int count[3], count_max[3];

        count[0] = -nth_omp_min;
        count[1] =  nth_omp_max;
        count[2] =  willUsePhysicalGpu;

        MPI_Allreduce(count, count_max, 3, MPI_INT, MPI_MAX, cr->mpi_comm_mysim);

        /* In case of an inhomogeneous run setup we use the maximum counts */
        nth_omp_min        = -count_max[0];
        nth_omp_max        =  count_max[1];
        anyRankIsUsingGpus = count_max[2] > 0;
    }

    int nthreads_omp_mpi_ok_min;

    if (!anyRankIsUsingGpus)
    {
        nthreads_omp_mpi_ok_min = nthreads_omp_mpi_ok_min_cpu;
    }
    else
    {
        /* With GPUs we set the minimum number of OpenMP threads to 2 to catch
         * cases where the user specifies #ranks == #cores.
         */
        nthreads_omp_mpi_ok_min = nthreads_omp_mpi_ok_min_gpu;
    }

    if (DOMAINDECOMP(cr) && cr->nnodes > 1)
    {
        if (nth_omp_max < nthreads_omp_mpi_ok_min ||
            nth_omp_max > nthreads_omp_mpi_ok_max)
        {
            /* Note that we print target_max here, not ok_max */
            sprintf(buf, "Your choice of number of MPI ranks and amount of resources results in using %d OpenMP threads per rank, which is most likely inefficient. The optimum is usually between %d and %d threads per rank.",
                    nth_omp_max,
                    nthreads_omp_mpi_ok_min,
                    nthreads_omp_mpi_target_max);

            if (bNtOmpOptionSet)
            {
                GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted("NOTE: %s", buf);
            }
            else
            {
                /* This fatal error, and the one below, is nasty, but it's
                 * probably the only way to ensure that all users don't waste
                 * a lot of resources, since many users don't read logs/stderr.
                 */
                gmx_fatal(FARGS, "%s If you want to run with this setup, specify the -ntomp option. But we suggest to change the number of MPI ranks%s.", buf, mpi_option);
            }
        }
    }
    else
    {
        const gmx::CpuInfo &cpuInfo = *hwinfo->cpuInfo;

        /* No domain decomposition (or only one domain) */
        if (nth_omp_max > nthreads_omp_faster(cpuInfo, anyRankIsUsingGpus))
        {
            /* To arrive here, the user/system set #ranks and/or #OMPthreads */
            gmx_bool bEnvSet;
            char     buf2[256];

            bEnvSet = (getenv("OMP_NUM_THREADS") != nullptr);

            if (bNtOmpOptionSet || bEnvSet)
            {
                sprintf(buf2, "You requested %d OpenMP threads", nth_omp_max);
            }
            else
            {
                sprintf(buf2, "Your choice of %d MPI rank%s and the use of %d total threads %sleads to the use of %d OpenMP threads",
                        cr->nnodes + cr->npmenodes,
                        cr->nnodes + cr->npmenodes == 1 ? "" : "s",
                        numTotalThreads > 0 ? numTotalThreads : hwinfo->nthreads_hw_avail,
                        hwinfo->nphysicalnode > 1 ? "on a node " : "",
                        nth_omp_max);
            }
            sprintf(buf, "%s, whereas we expect the optimum to be with more MPI ranks with %d to %d OpenMP threads.",
                    buf2, nthreads_omp_mpi_ok_min, nthreads_omp_mpi_target_max);

            /* We can not quit with a fatal error when OMP_NUM_THREADS is set
             * with different values per rank or node, since in that case
             * the user can not set -ntomp to override the error.
             */
            if (bNtOmpOptionSet || (bEnvSet && nth_omp_min != nth_omp_max))
            {
                GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted("NOTE: %s", buf);
            }
            else
            {
                gmx_fatal(FARGS, "%s If you want to run with this many OpenMP threads, specify the -ntomp option. But we suggest to increase the number of MPI ranks%s.", buf, mpi_option);
            }
        }
    }
#else /* GMX_OPENMP && GMX_MPI */
      /* No OpenMP and/or MPI: it doesn't make much sense to check */
    GMX_UNUSED_VALUE(bNtOmpOptionSet);
    GMX_UNUSED_VALUE(numTotalThreads);
    GMX_UNUSED_VALUE(willUsePhysicalGpu);
    GMX_UNUSED_VALUE(cr);
    /* Check if we have more than 1 physical core, if detected,
     * or more than 1 hardware thread if physical cores were not detected.
     */
    if (!GMX_OPENMP && !GMX_MPI && hwinfo->hardwareTopology->numberOfCores() > 1)
    {
        GMX_LOG(mdlog.warning).asParagraph().appendText("NOTE: GROMACS was compiled without OpenMP and (thread-)MPI support, can only use a single CPU core");
    }
#endif /* GMX_OPENMP && GMX_MPI */
}


//! Dump a \c hw_opt to \c fp.
static void print_hw_opt(FILE *fp, const gmx_hw_opt_t *hw_opt)
{
    fprintf(fp, "hw_opt: nt %d ntmpi %d ntomp %d ntomp_pme %d gpu_id '%s'\n",
            hw_opt->nthreads_tot,
            hw_opt->nthreads_tmpi,
            hw_opt->nthreads_omp,
            hw_opt->nthreads_omp_pme,
            hw_opt->gpuIdTaskAssignment.c_str());
}

void check_and_update_hw_opt_1(gmx_hw_opt_t    *hw_opt,
                               const t_commrec *cr,
                               int              nPmeRanks)
{
    /* Currently hw_opt only contains default settings or settings supplied
     * by the user on the command line.
     */
    if (hw_opt->nthreads_omp < 0)
    {
        gmx_fatal(FARGS, "The number of OpenMP threads supplied on the command line is %d, which is negative and not allowed", hw_opt->nthreads_omp);
    }

    /* Check for OpenMP settings stored in environment variables, which can
     * potentially be different on different MPI ranks.
     */
    gmx_omp_nthreads_read_env(&hw_opt->nthreads_omp, SIMMASTER(cr));

    /* Check restrictions on the user supplied options before modifying them.
     * TODO: Put the user values in a const struct and preserve them.
     */
#if !GMX_THREAD_MPI
    if (hw_opt->nthreads_tot > 0)
    {
        gmx_fatal(FARGS, "Setting the total number of threads is only supported with thread-MPI and GROMACS was compiled without thread-MPI");
    }
    if (hw_opt->nthreads_tmpi > 0)
    {
        gmx_fatal(FARGS, "Setting the number of thread-MPI ranks is only supported with thread-MPI and GROMACS was compiled without thread-MPI");
    }
#endif

    if (bHasOmpSupport)
    {
        /* Check restrictions on PME thread related options set by the user */

        if (hw_opt->nthreads_omp_pme > 0 && hw_opt->nthreads_omp <= 0)
        {
            gmx_fatal(FARGS, "You need to specify -ntomp in addition to -ntomp_pme");
        }

        if (hw_opt->nthreads_omp_pme >= 1 &&
            hw_opt->nthreads_omp_pme != hw_opt->nthreads_omp &&
            nPmeRanks <= 0)
        {
            /* This can result in a fatal error on many MPI ranks,
             * but since the thread count can differ per rank,
             * we can't easily avoid this.
             */
            gmx_fatal(FARGS, "You need to explicitly specify the number of PME ranks (-npme) when using different number of OpenMP threads for PP and PME ranks");
        }
    }
    else
    {
        /* GROMACS was configured without OpenMP support */

        if (hw_opt->nthreads_omp > 1 || hw_opt->nthreads_omp_pme > 1)
        {
            gmx_fatal(FARGS, "More than 1 OpenMP thread requested, but GROMACS was compiled without OpenMP support");
        }
        hw_opt->nthreads_omp     = 1;
        hw_opt->nthreads_omp_pme = 1;
    }

    if (hw_opt->nthreads_tot > 0 && hw_opt->nthreads_omp_pme <= 0)
    {
        /* We have the same number of OpenMP threads for PP and PME ranks,
         * thus we can perform several consistency checks.
         */
        if (hw_opt->nthreads_tmpi > 0 &&
            hw_opt->nthreads_omp > 0 &&
            hw_opt->nthreads_tot != hw_opt->nthreads_tmpi*hw_opt->nthreads_omp)
        {
            gmx_fatal(FARGS, "The total number of threads requested (%d) does not match the thread-MPI ranks (%d) times the OpenMP threads (%d) requested",
                      hw_opt->nthreads_tot, hw_opt->nthreads_tmpi, hw_opt->nthreads_omp);
        }

        if (hw_opt->nthreads_tmpi > 0 &&
            hw_opt->nthreads_tot % hw_opt->nthreads_tmpi != 0)
        {
            gmx_fatal(FARGS, "The total number of threads requested (%d) is not divisible by the number of thread-MPI ranks requested (%d)",
                      hw_opt->nthreads_tot, hw_opt->nthreads_tmpi);
        }

        if (hw_opt->nthreads_omp > 0 &&
            hw_opt->nthreads_tot % hw_opt->nthreads_omp != 0)
        {
            gmx_fatal(FARGS, "The total number of threads requested (%d) is not divisible by the number of OpenMP threads requested (%d)",
                      hw_opt->nthreads_tot, hw_opt->nthreads_omp);
        }
    }

    if (hw_opt->nthreads_tot > 0)
    {
        if (hw_opt->nthreads_omp > hw_opt->nthreads_tot)
        {
            gmx_fatal(FARGS, "You requested %d OpenMP threads with %d total threads. Choose a total number of threads that is a multiple of the number of OpenMP threads.",
                      hw_opt->nthreads_omp, hw_opt->nthreads_tot);
        }

        if (hw_opt->nthreads_tmpi > hw_opt->nthreads_tot)
        {
            gmx_fatal(FARGS, "You requested %d thread-MPI ranks with %d total threads. Choose a total number of threads that is a multiple of the number of thread-MPI ranks.",
                      hw_opt->nthreads_tmpi, hw_opt->nthreads_tot);
        }
    }

    if (debug)
    {
        print_hw_opt(debug, hw_opt);
    }

    /* Asserting this simplifies the hardware resource division later
     * on. */
    GMX_RELEASE_ASSERT(!(hw_opt->nthreads_omp_pme >= 1 && hw_opt->nthreads_omp <= 0),
                       "PME thread count should only be set when the normal thread count is also set");
}

void check_and_update_hw_opt_2(gmx_hw_opt_t *hw_opt,
                               int           cutoff_scheme)
{
    if (cutoff_scheme == ecutsGROUP)
    {
        /* We only have OpenMP support for PME only nodes */
        if (hw_opt->nthreads_omp > 1)
        {
            gmx_fatal(FARGS, "OpenMP threads have been requested with cut-off scheme %s, but these are only supported with cut-off scheme %s",
                      ecutscheme_names[cutoff_scheme],
                      ecutscheme_names[ecutsVERLET]);
        }
        hw_opt->nthreads_omp = 1;
    }
}

void check_and_update_hw_opt_3(gmx_hw_opt_t *hw_opt)
{
#if GMX_THREAD_MPI
    GMX_RELEASE_ASSERT(hw_opt->nthreads_tmpi >= 1, "Must have at least one thread-MPI rank");

    /* If the user set the total number of threads on the command line
     * and did not specify the number of OpenMP threads, set the latter here.
     */
    if (hw_opt->nthreads_tot > 0 && hw_opt->nthreads_omp <= 0)
    {
        hw_opt->nthreads_omp = hw_opt->nthreads_tot/hw_opt->nthreads_tmpi;

        if (!bHasOmpSupport && hw_opt->nthreads_omp > 1)
        {
            gmx_fatal(FARGS, "You (indirectly) asked for OpenMP threads by setting -nt > -ntmpi, but GROMACS was compiled without OpenMP support");
        }
    }
#endif

    GMX_RELEASE_ASSERT(bHasOmpSupport || hw_opt->nthreads_omp == 1, "Without OpenMP support, only one thread per rank can be used");

    /* We are done with updating nthreads_omp, we can set nthreads_omp_pme */
    if (hw_opt->nthreads_omp_pme <= 0 && hw_opt->nthreads_omp > 0)
    {
        hw_opt->nthreads_omp_pme = hw_opt->nthreads_omp;
    }

    if (debug)
    {
        print_hw_opt(debug, hw_opt);
    }
}
