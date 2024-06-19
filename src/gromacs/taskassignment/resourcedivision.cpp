/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
/*! \internal \file
 * \brief Defines utility functionality for dividing resources and
 * checking for consistency and usefulness.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */

#include "gmxpre.h"

#include "gromacs/taskassignment/resourcedivision.h"

#include "config.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/ewald/pme.h"
#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"
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

/*! \brief The minimum number of atoms per thread-MPI thread when GPUs
 * are present. With fewer atoms than this, the number of thread-MPI
 * ranks will get lowered.
 */
static constexpr int min_atoms_per_mpi_thread = 90;
/*! \brief The minimum number of atoms per GPU with thread-MPI
 * active. With fewer atoms than this, the number of thread-MPI ranks
 * will get lowered.
 */
static constexpr int min_atoms_per_gpu = 900;

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
constexpr int nthreads_omp_faster_default   = 8;
constexpr int nthreads_omp_faster_Nehalem   = 12;
constexpr int nthreads_omp_faster_Intel_AVX = 16;
constexpr int nthreads_omp_faster_AMD_Ryzen = 16;
/* For CPU only runs the fastest options are usually MPI or OpenMP only.
 * With one GPU, using MPI only is almost never optimal, so we need to
 * compare running pure OpenMP with combined MPI+OpenMP. This means higher
 * OpenMP threads counts can still be ok. Multiplying the numbers above
 * by a factor of 2 seems to be a good estimate.
 */
constexpr int nthreads_omp_faster_gpu_fac = 2;

/* This is the case with MPI (2 or more MPI PP ranks).
 * By default we will terminate with a fatal error when more than 8
 * OpenMP thread are (indirectly) requested, since using less threads
 * nearly always results in better performance.
 * With thread-mpi and multiple GPUs or one GPU and too many threads
 * we first try 6 OpenMP threads and then less until the number of MPI ranks
 * is divisible by the number of GPUs.
 */
constexpr int nthreads_omp_mpi_ok_max     = 8;
constexpr int nthreads_omp_mpi_ok_min_cpu = 1;
constexpr int nthreads_omp_mpi_ok_min_gpu = 2;
constexpr int nthreads_omp_mpi_target_max = 6;

// Too many ranks per GPU can lead to large overhead so we cap the
// tMPI rank count choosen automatically.
constexpr int c_maxAutoTmpiRanksPerGpu = 4;
/**@}*/

/*! \brief Returns the maximum OpenMP thread count for which using a single MPI rank
 * should be faster than using multiple ranks with the same total thread count.
 */
static int nthreads_omp_faster(const gmx::CpuInfo& cpuInfo, gmx_bool bUseGPU)
{
    int nth;

    if (cpuInfo.vendor() == gmx::CpuInfo::Vendor::Intel && cpuInfo.feature(gmx::CpuInfo::Feature::X86_Avx))
    {
        nth = nthreads_omp_faster_Intel_AVX;
    }
    else if (gmx::cpuIsX86Nehalem(cpuInfo))
    {
        // Intel Nehalem
        nth = nthreads_omp_faster_Nehalem;
    }
    else if ((cpuInfo.vendor() == gmx::CpuInfo::Vendor::Amd && cpuInfo.family() >= 23)
             || cpuInfo.vendor() == gmx::CpuInfo::Vendor::Hygon)
    {
        // AMD Ryzen || Hygon Dhyana
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
gmx_unused static int nthreads_omp_efficient_max(int gmx_unused nrank, const gmx::CpuInfo& cpuInfo, gmx_bool bUseGPU)
{
    if (GMX_OPENMP && GMX_MPI && (nrank > 1))
    {
        return nthreads_omp_mpi_ok_max;
    }
    else
    {
        return nthreads_omp_faster(cpuInfo, bUseGPU);
    }
}

/*! \brief Return the number of thread-MPI ranks to use.
 * This is chosen such that we can always obey our own efficiency checks.
 */
gmx_unused static int get_tmpi_omp_thread_division(const gmx_hw_info_t* hwinfo,
                                                   const gmx_hw_opt_t&  hw_opt,
                                                   int                  nthreads_tot,
                                                   int                  ngpu)
{
    int                 nrank;
    const gmx::CpuInfo& cpuInfo = *hwinfo->cpuInfo;

    GMX_RELEASE_ASSERT(nthreads_tot > 0, "There must be at least one thread per rank");

    /* There are no separate PME nodes here, as we ensured in
     * check_and_update_hw_opt that nthreads_tmpi>0 with PME nodes
     * and a conditional ensures we would not have ended up here.
     * Note that separate PME nodes might be switched on later.
     */
    if (ngpu > 0)
    {
        if (hw_opt.nthreads_omp > 0)
        {
            /* In this case it is unclear if we should use 1 rank per GPU
             * or more or less, so we require also setting the number of ranks.
             */
            gmx_fatal(FARGS,
                      "When using GPUs, setting the number of OpenMP threads without specifying "
                      "the number "
                      "of ranks can lead to conflicting demands. Please specify the number of "
                      "thread-MPI ranks "
                      "as well (option -ntmpi).");
        }

        nrank = ngpu;

        /* When the user sets nthreads_omp, we can end up oversubscribing CPU cores
         * if we simply start as many ranks as GPUs. To avoid this, we start as few
         * tMPI ranks as necessary to avoid oversubscription and instead leave GPUs idle.
         * If the user does not set the number of OpenMP threads, nthreads_omp==0 and
         * this code has no effect.
         */
        GMX_RELEASE_ASSERT(hw_opt.nthreads_omp >= 0,
                           "nthreads_omp is negative, but previous checks should "
                           "have prevented this");
        while (nrank * hw_opt.nthreads_omp > hwinfo->hardwareTopology->maxThreads() && nrank > 1)
        {
            nrank--;
        }

        if (nthreads_tot < nrank)
        {
            /* #thread < #gpu is very unlikely, but if so: waste gpu(s) */
            nrank = nthreads_tot;
        }
        else if (nthreads_tot > nthreads_omp_faster(cpuInfo, ngpu > 0)
                 || (ngpu > 1 && nthreads_tot / ngpu > nthreads_omp_mpi_target_max))
        {
            /* The high OpenMP thread count will likely result in sub-optimal
             * performance. Increase the rank count to reduce the thread count
             * per rank. This will lead to GPU sharing by MPI ranks/threads.
             */
            int nshare;

            // Increase the rank count as long as have we more than nthreads_omp_mpi_target_max OpenMP
            // threads per rank or either the number of hardware threads is not
            // divisible by the rank count or we would have >4 ranks per GPU (which is inefficient).
            // Don't go below nthreads_omp_mpi_ok_min_gpu OpenMP threads/rank.
            nshare = 1;
            do
            {
                nshare++;
                nrank = ngpu * nshare;
            } while ((nthreads_tot / nrank > nthreads_omp_mpi_target_max && nshare < c_maxAutoTmpiRanksPerGpu)
                     || (nthreads_tot / (ngpu * (nshare + 1)) >= nthreads_omp_mpi_ok_min_gpu
                         && nthreads_tot % nrank != 0));
        }
    }
    else if (hw_opt.nthreads_omp > 0)
    {
        /* Here we could oversubscribe, when we do, we issue a warning later */
        nrank = std::max(1, nthreads_tot / hw_opt.nthreads_omp);
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

//! Return whether hyper threading is used on ALL cores.
static bool gmxSmtIsUsedOnAllCores(const gmx::HardwareTopology& hwTop)
{
    if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Basic)
    {
        std::size_t minSmt = 999999999;
        std::size_t maxSmt = 0;

        for (const auto& p : hwTop.machine().packages)
        {
            for (const auto& c : p.cores)
            {
                minSmt = std::min(minSmt, c.processingUnits.size());
                maxSmt = std::max(maxSmt, c.processingUnits.size());
            }
        }
        if (minSmt == maxSmt && minSmt > 1)
        {
            return true;
        }
    }
    return false;
}

namespace
{

//! Handles checks for algorithms that must use a single rank.
class SingleRankChecker
{
public:
    SingleRankChecker() : value_(false) {}
    /*! \brief Call this function for each possible condition
        under which a single rank is required, along with a string
        describing the constraint when it is applied. */
    void applyConstraint(bool condition, const char* description)
    {
        if (condition)
        {
            value_ = true;
            reasons_.push_back(gmx::formatString("%s only supports a single rank.", description));
        }
    }
    //! After applying any conditions, is a single rank required?
    bool mustUseOneRank() const { return value_; }
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
int get_nthreads_mpi(const gmx_hw_info_t* hwinfo,
                     gmx_hw_opt_t*        hw_opt,
                     const int            numDevicesToUse,
                     bool                 nonbondedOnGpu,
                     bool                 pmeOnGpu,
                     const t_inputrec*    inputrec,
                     const gmx_mtop_t&    mtop,
                     const gmx::MDLogger& mdlog,
                     bool                 doMembed)
{
    int nthreads_hw, nthreads_tot_max, nrank, ngpu;
    int min_atoms_per_mpi_rank;

    const gmx::CpuInfo&          cpuInfo = *hwinfo->cpuInfo;
    const gmx::HardwareTopology& hwTop   = *hwinfo->hardwareTopology;

    if (pmeOnGpu)
    {
        GMX_RELEASE_ASSERT((usingPme(inputrec->coulombtype) || usingLJPme(inputrec->vdwtype))
                                   && pme_gpu_supports_build(nullptr)
                                   && pme_gpu_supports_input(*inputrec, nullptr),
                           "PME can't be on GPUs unless we are using PME");

        // PME on GPUs supports a single PME rank with PP running on the same or few other ranks.
        // For now, let's treat separate PME GPU rank as opt-in.
        if (hw_opt->nthreads_tmpi < 1)
        {
            return 1;
        }
    }

    {
        /* Check if an algorithm does not support parallel simulation.  */
        // TODO This might work better if e.g. implemented algorithms
        // had to define a function that returns such requirements,
        // and a description string.
        SingleRankChecker checker;
        checker.applyConstraint(inputrec->eI == IntegrationAlgorithm::LBFGS, "L-BFGS minimization");
        checker.applyConstraint(inputrec->coulombtype == CoulombInteractionType::Ewald,
                                "Plain Ewald electrostatics");
        checker.applyConstraint(doMembed, "Membrane embedding");
        bool useOrientationRestraints = (gmx_mtop_ftype_count(mtop, F_ORIRES) > 0);
        checker.applyConstraint(useOrientationRestraints, "Orientation restraints");
        if (checker.mustUseOneRank())
        {
            std::string message = checker.getMessage();
            if (hw_opt->nthreads_tmpi > 1)
            {
                gmx_fatal(FARGS,
                          "%s However, you asked for more than 1 thread-MPI rank, so mdrun cannot "
                          "continue. "
                          "Choose a single rank, or a different algorithm.",
                          message.c_str());
            }
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted("%s Choosing to use only a single thread-MPI rank.",
                                         message.c_str());
            return 1;
        }
    }

    if (hw_opt->nthreads_tmpi > 0)
    {
        /* Trivial, return the user's choice right away */
        return hw_opt->nthreads_tmpi;
    }

    // Now implement automatic selection of number of thread-MPI ranks
    nthreads_hw = hwinfo->hardwareTopology->maxThreads();

    if (nthreads_hw <= 0)
    {
        /* This should normally not happen, but if it does, we handle it */
        gmx_fatal(FARGS,
                  "The number of available hardware threads can not be detected, please specify "
                  "the number of "
                  "MPI ranks and the number of OpenMP threads (if supported) manually with options "
                  "-ntmpi and -ntomp, respectively");
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

    /* nonbondedOnGpu might be false e.g. because this simulation
     * is a rerun with energy groups. */
    ngpu = (nonbondedOnGpu ? numDevicesToUse : 0);

    nrank = get_tmpi_omp_thread_division(hwinfo, *hw_opt, nthreads_tot_max, ngpu);

    if (inputrec->eI == IntegrationAlgorithm::NM || EI_TPI(inputrec->eI))
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

    if (mtop.natoms / nrank < min_atoms_per_mpi_rank)
    {
        int nrank_new;

        /* the rank number was chosen automatically, but there are too few
           atoms per rank, so we need to reduce the rank count */
        nrank_new = std::max(1, mtop.natoms / min_atoms_per_mpi_rank);

        /* Avoid partial use of Hyper-Threading */
        if (gmxSmtIsUsedOnAllCores(hwTop) && nrank_new > nthreads_hw / 2 && nrank_new < nthreads_hw)
        {
            nrank_new = nthreads_hw / 2;
        }

        /* If the user specified the total thread count, ensure this is
         * divisible by the number of ranks.
         * It is quite likely that we have too many total threads compared
         * to the size of the system, but if the user asked for this many
         * threads we should respect that.
         */
        while (hw_opt->nthreads_tot > 0 && hw_opt->nthreads_tot % nrank_new != 0)
        {
            nrank_new--;
        }

        /* Avoid large prime numbers in the rank count */
        if (nrank_new >= 6)
        {
            /* Use only 6,8,10 with additional factors of 2 */
            int fac;

            fac = 2;
            while (3 * fac * 2 <= nrank_new)
            {
                fac *= 2;
            }

            nrank_new = (nrank_new / fac) * fac;
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

        if (ngpu > 0 && (nrank_new % ngpu) != 0)
        {
            /* If we use GPUs, the number of ranks must be divisible by the number of GPUs,
             * unless the GPUs are very different (and if they are, user should manually
             * select the parallelization scheme).
             * Rounding down the number of ranks and limiting the total count per GPU.
             * */
            if (nrank_new > ngpu)
            {
                nrank_new = std::min(nrank_new / ngpu, c_maxAutoTmpiRanksPerGpu) * ngpu;
            }
            else
            {
                nrank_new = ngpu;
            }
        }

        nrank = nrank_new;

        /* We reduced the number of tMPI ranks, which means we might violate
         * our own efficiency checks if we simply use all hardware threads.
         */
        if (GMX_OPENMP && hw_opt->nthreads_omp <= 0 && hw_opt->nthreads_tot <= 0)
        {
            /* The user set neither the total nor the OpenMP thread count,
             * we should use all hardware threads, unless we will violate
             * our own efficiency limitation on the thread count.
             */
            int nt_omp_max;

            nt_omp_max = nthreads_omp_efficient_max(nrank, cpuInfo, ngpu >= 1);

            if (nrank * nt_omp_max < hwinfo->hardwareTopology->maxThreads())
            {
                /* Limit the number of OpenMP threads to start */
                hw_opt->nthreads_omp = nt_omp_max;
            }
        }

        fprintf(stderr, "\n");
        fprintf(stderr, "NOTE: Parallelization is limited by the small number of atoms,\n");
        fprintf(stderr, "      only starting %d thread-MPI ranks.\n", nrank);
        fprintf(stderr,
                "      You can use the -nt and/or -ntmpi option to optimize the number of "
                "threads.\n\n");
    }

    return nrank;
}


void check_resource_division_efficiency(const gmx_hw_info_t* hwinfo,
                                        bool                 willUsePhysicalGpu,
                                        t_commrec*           cr,
                                        const gmx::MDLogger& mdlog)
{
    GMX_UNUSED_VALUE(hwinfo);
#if GMX_OPENMP && GMX_MPI

    /* This function should be called after thread-MPI (when configured) and
     * OpenMP have been initialized. Check that here.
     */
    if (GMX_THREAD_MPI)
    {
        GMX_RELEASE_ASSERT(nthreads_omp_faster_default >= nthreads_omp_mpi_ok_max,
                           "Inconsistent OpenMP thread count default values");
    }
    GMX_RELEASE_ASSERT(gmx_omp_nthreads_get(ModuleMultiThread::Default) >= 1,
                       "Must have at least one OpenMP thread");

    int nth_omp_max = gmx_omp_nthreads_get(ModuleMultiThread::Default);

    bool anyRankIsUsingGpus = willUsePhysicalGpu;
    /* Thread-MPI seems to have a bug with reduce on 1 node, so use a cond. */
    if (cr->nnodes > 1)
    {
        std::array<int, 2> count, count_max;

        count[0] = nth_omp_max;
        count[1] = int(willUsePhysicalGpu);

        MPI_Allreduce(count.data(), count_max.data(), count.size(), MPI_INT, MPI_MAX, cr->mpi_comm_mysim);

        /* In case of an inhomogeneous run setup we use the maximum counts */
        nth_omp_max        = count_max[0];
        anyRankIsUsingGpus = count_max[1] > 0;
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

    if (cr && (cr->nnodes > 1) && !anyRankIsUsingGpus)
    {
        if (nth_omp_max < nthreads_omp_mpi_ok_min || nth_omp_max > nthreads_omp_mpi_ok_max)
        {
            auto msg = gmx::formatString(
                    "Note: Your choice of number of MPI ranks and amount of resources results in "
                    "using %d OpenMP threads per rank, which is most likely inefficient. "
                    "The optimum is usually between %d and %d threads per rank.",
                    nth_omp_max,
                    nthreads_omp_mpi_ok_min,
                    nthreads_omp_mpi_ok_max);

            GMX_LOG(mdlog.warning).asParagraph().appendText(msg);
        }
    }
#else  // !GMX_OPENMP || ! GMX_MPI
    GMX_UNUSED_VALUE(willUsePhysicalGpu);
    GMX_UNUSED_VALUE(cr);
    GMX_UNUSED_VALUE(nthreads_omp_mpi_ok_max);
    GMX_UNUSED_VALUE(nthreads_omp_mpi_ok_min_cpu);

    // Since even an Apple watch comes with 2 cores today, we can probably safely assume we
    // aren't running on a single-core device and warn the user if we can't use many threads.
    if (!GMX_OPENMP && !GMX_MPI)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "NOTE: GROMACS was compiled without OpenMP and (thread-)MPI support, can "
                        "only use a single CPU core");
    }
#endif // end GMX_OPENMP && GMX_MPI
}


//! Dump a \c hw_opt to \c fp.
static void print_hw_opt(FILE* fp, const gmx_hw_opt_t* hw_opt)
{
    fprintf(fp,
            "hw_opt: nt %d ntmpi %d ntomp %d ntomp_pme %d gpu_id '%s' gputasks '%s'\n",
            hw_opt->nthreads_tot,
            hw_opt->nthreads_tmpi,
            hw_opt->nthreads_omp,
            hw_opt->nthreads_omp_pme,
            hw_opt->devicesSelectedByUser.c_str(),
            hw_opt->userGpuTaskAssignment.c_str());
}

void checkAndUpdateHardwareOptions(const gmx::MDLogger& mdlog,
                                   gmx_hw_opt_t*        hw_opt,
                                   const bool           isSimulationMainRank,
                                   const int            nPmeRanks,
                                   const t_inputrec*    inputrec)
{
    /* Currently hw_opt only contains default settings or settings supplied
     * by the user on the command line.
     */
    if (hw_opt->nthreads_omp < 0)
    {
        gmx_fatal(FARGS,
                  "The number of OpenMP threads supplied on the command line is %d, which is "
                  "negative "
                  "and not allowed",
                  hw_opt->nthreads_omp);
    }

    /* Check for OpenMP settings stored in environment variables, which can
     * potentially be different on different MPI ranks.
     */
    gmx_omp_nthreads_read_env(mdlog, &hw_opt->nthreads_omp);

    /* Check restrictions on the user supplied options before modifying them.
     * TODO: Put the user values in a const struct and preserve them.
     */
    if (!GMX_THREAD_MPI)
    {

        if (hw_opt->nthreads_tot > 0)
        {
            gmx_fatal(FARGS,
                      "Setting the total number of threads is only supported with thread-MPI and "
                      "GROMACS was "
                      "compiled without thread-MPI");
        }
        if (hw_opt->nthreads_tmpi > 0)
        {
            gmx_fatal(FARGS,
                      "Setting the number of thread-MPI ranks is only supported with thread-MPI "
                      "and GROMACS was "
                      "compiled without thread-MPI");
        }
    }

    /* With thread-MPI we need to handle TPI and #OpenMP-threads=auto early,
     * so we can parallelize using MPI only. The general check is done later.
     */
    if (GMX_THREAD_MPI && isSimulationMainRank)
    {
        GMX_RELEASE_ASSERT(inputrec, "Expect a valid inputrec");
        if (EI_TPI(inputrec->eI) && hw_opt->nthreads_omp == 0)
        {
            hw_opt->nthreads_omp = 1;
        }
    }
    /* With thread-MPI the main thread sets hw_opt->totNumThreadsIsAuto.
     * The other threads receive a partially processed hw_opt from the main
     * thread and should not set hw_opt->totNumThreadsIsAuto again.
     */
    if (!GMX_THREAD_MPI || isSimulationMainRank)
    {
        /* Check if mdrun is free to choose the total number of threads */
        hw_opt->totNumThreadsIsAuto = (hw_opt->nthreads_omp == 0 && hw_opt->nthreads_omp_pme == 0
                                       && hw_opt->nthreads_tot == 0);
    }

    if (GMX_OPENMP)
    {
        /* Check restrictions on PME thread related options set by the user */

        if (hw_opt->nthreads_omp_pme > 0 && hw_opt->nthreads_omp <= 0)
        {
            gmx_fatal(FARGS, "You need to specify -ntomp in addition to -ntomp_pme");
        }

        if (hw_opt->nthreads_omp_pme >= 1 && hw_opt->nthreads_omp_pme != hw_opt->nthreads_omp
            && nPmeRanks <= 0)
        {
            /* This can result in a fatal error on many MPI ranks,
             * but since the thread count can differ per rank,
             * we can't easily avoid this.
             */
            gmx_fatal(FARGS,
                      "You need to explicitly specify the number of PME ranks (-npme) when using "
                      "different numbers of OpenMP threads for PP and PME ranks");
        }
    }
    else
    {
        /* GROMACS was configured without OpenMP support */

        if (hw_opt->nthreads_omp > 1 || hw_opt->nthreads_omp_pme > 1)
        {
            gmx_fatal(FARGS,
                      "More than 1 OpenMP thread requested, but GROMACS was compiled without "
                      "OpenMP support");
        }
        hw_opt->nthreads_omp     = 1;
        hw_opt->nthreads_omp_pme = 1;
    }

    if (hw_opt->nthreads_tot > 0 && hw_opt->nthreads_omp_pme <= 0)
    {
        /* We have the same number of OpenMP threads for PP and PME ranks,
         * thus we can perform several consistency checks.
         */
        if (hw_opt->nthreads_tmpi > 0 && hw_opt->nthreads_omp > 0
            && hw_opt->nthreads_tot != hw_opt->nthreads_tmpi * hw_opt->nthreads_omp)
        {
            gmx_fatal(FARGS,
                      "The total number of threads requested (%d) does not match the thread-MPI "
                      "ranks (%d) "
                      "times the OpenMP threads (%d) requested",
                      hw_opt->nthreads_tot,
                      hw_opt->nthreads_tmpi,
                      hw_opt->nthreads_omp);
        }

        if (hw_opt->nthreads_tmpi > 0 && hw_opt->nthreads_tot % hw_opt->nthreads_tmpi != 0)
        {
            gmx_fatal(FARGS,
                      "The total number of threads requested (%d) is not divisible by the number "
                      "of thread-MPI "
                      "ranks requested (%d)",
                      hw_opt->nthreads_tot,
                      hw_opt->nthreads_tmpi);
        }

        if (hw_opt->nthreads_omp > 0 && hw_opt->nthreads_tot % hw_opt->nthreads_omp != 0)
        {
            gmx_fatal(FARGS,
                      "The total number of threads requested (%d) is not divisible by the number "
                      "of OpenMP "
                      "threads requested (%d)",
                      hw_opt->nthreads_tot,
                      hw_opt->nthreads_omp);
        }
    }

    if (hw_opt->nthreads_tot > 0)
    {
        if (hw_opt->nthreads_omp > hw_opt->nthreads_tot)
        {
            gmx_fatal(FARGS,
                      "You requested %d OpenMP threads with %d total threads. Choose a total "
                      "number of threads "
                      "that is a multiple of the number of OpenMP threads.",
                      hw_opt->nthreads_omp,
                      hw_opt->nthreads_tot);
        }

        if (hw_opt->nthreads_tmpi > hw_opt->nthreads_tot)
        {
            gmx_fatal(FARGS,
                      "You requested %d thread-MPI ranks with %d total threads. Choose a total "
                      "number of "
                      "threads that is a multiple of the number of thread-MPI ranks.",
                      hw_opt->nthreads_tmpi,
                      hw_opt->nthreads_tot);
        }
    }

    if (GMX_THREAD_MPI && nPmeRanks > 0 && hw_opt->nthreads_tmpi <= 0)
    {
        gmx_fatal(FARGS,
                  "You need to explicitly specify the number of MPI threads (-ntmpi) when using "
                  "separate PME ranks");
    }

    if (debug)
    {
        print_hw_opt(debug, hw_opt);
    }

    /* Asserting this simplifies the hardware resource division later
     * on. */
    GMX_RELEASE_ASSERT(
            !(hw_opt->nthreads_omp_pme >= 1 && hw_opt->nthreads_omp <= 0),
            "PME thread count should only be set when the normal thread count is also set");
}

void checkAndUpdateRequestedNumOpenmpThreads(gmx_hw_opt_t*         hw_opt,
                                             const gmx_hw_info_t&  hwinfo,
                                             const t_commrec*      cr,
                                             const gmx_multisim_t* ms,
                                             int                   numRanksOnThisNode,
                                             PmeRunMode            pmeRunMode,
                                             const gmx_mtop_t&     mtop,
                                             const t_inputrec&     inputrec)
{
    if (EI_TPI(inputrec.eI))
    {
        if (hw_opt->nthreads_omp > 1)
        {
            gmx_fatal(FARGS,
                      "You requested OpenMP parallelization, which is not supported with TPI.");
        }
        hw_opt->nthreads_omp = 1;
    }

    if (GMX_THREAD_MPI)
    {

        GMX_RELEASE_ASSERT(hw_opt->nthreads_tmpi >= 1, "Must have at least one thread-MPI rank");

        /* If the user set the total number of threads on the command line
         * and did not specify the number of OpenMP threads, set the latter here.
         */
        if (hw_opt->nthreads_tot > 0 && hw_opt->nthreads_omp <= 0)
        {
            hw_opt->nthreads_omp = hw_opt->nthreads_tot / hw_opt->nthreads_tmpi;

            if (!GMX_OPENMP && hw_opt->nthreads_omp > 1)
            {
                gmx_fatal(FARGS,
                          "You (indirectly) asked for OpenMP threads by setting -nt > -ntmpi, but "
                          "GROMACS was "
                          "compiled without OpenMP support");
            }
        }
    }
    /* With both non-bonded and PME on GPU, the work left on the CPU is often
     * (much) slower with SMT than without SMT. This is mostly the case with
     * few atoms per core. Thus, if the number of threads is set to auto,
     * we turn off SMT in that case. Note that PME on GPU implies that also
     * the non-bonded are computed on the GPU.
     * We only need to do this when the number of hardware theads is larger
     * than the number of cores. Note that a queuing system could limit
     * the number of hardware threads available, but we are not trying to be
     * too smart here in that case.
     */
    /* The thread reduction and synchronization costs go up roughy quadratically
     * with the threads count, so we apply a threshold quadratic in #cores.
     * Also more cores per GPU usually means the CPU gets faster than the GPU.
     * The number 1000 atoms per core^2 is a reasonable threshold
     * for Intel x86 and AMD Threadripper.
     */
    constexpr int c_numAtomsPerCoreSquaredSmtThreshold = 1000;

    /* Prepare conditions for deciding if we should disable SMT.
     * We currently only limit SMT for simulations using a single rank.
     * TODO: Consider limiting also for multi-rank simulations.
     */
    bool canChooseNumOpenmpThreads      = (GMX_OPENMP && hw_opt->nthreads_omp <= 0);
    bool haveSmtSupport                 = gmxSmtIsUsedOnAllCores(*hwinfo.hardwareTopology);
    bool simRunsSingleRankNBAndPmeOnGpu = (cr->nnodes == 1 && pmeRunMode == PmeRunMode::GPU);

    if (canChooseNumOpenmpThreads && haveSmtSupport && simRunsSingleRankNBAndPmeOnGpu)
    {
        /* Note that the queing system might have limited us from using
         * all detected ncore_tot physical cores. We are currently not
         * checking for that here.
         */
        int numRanksTot     = cr->nnodes * (isMultiSim(ms) ? ms->numSimulations_ : 1);
        int numAtomsPerRank = mtop.natoms / cr->nnodes;
        int numCoresPerRank = hwinfo.ncore_tot / numRanksTot;
        if (numAtomsPerRank < c_numAtomsPerCoreSquaredSmtThreshold * gmx::square(numCoresPerRank))
        {
            // Don't run more than one thread per core
            int nCores = 0;
            for (const auto& p : hwinfo.hardwareTopology->machine().packages)
            {
                nCores += p.cores.size();
            }
            int maxThreads = std::min(nCores, hwinfo.hardwareTopology->maxThreads());

            hw_opt->nthreads_omp = std::max(1, maxThreads / numRanksOnThisNode);
        }
    }

    GMX_RELEASE_ASSERT(GMX_OPENMP || hw_opt->nthreads_omp == 1,
                       "Without OpenMP support, only one thread per rank can be used");

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

namespace gmx
{

void checkHardwareOversubscription(int                             numThreadsOnThisRank,
                                   int                             rank,
                                   const HardwareTopology&         hwTop,
                                   const PhysicalNodeCommunicator& comm,
                                   const MDLogger&                 mdlog)
{
    if (hwTop.supportLevel() < HardwareTopology::SupportLevel::LogicalProcessorCount)
    {
        /* There is nothing we can check */
        return;
    }

    int numRanksOnThisNode   = comm.size_;
    int numThreadsOnThisNode = numThreadsOnThisRank;
    /* Avoid MPI calls with uninitialized thread-MPI communicators */
    if (comm.size_ > 1)
    {
#if GMX_MPI
        /* Count the threads within this physical node */
        MPI_Allreduce(&numThreadsOnThisRank, &numThreadsOnThisNode, 1, MPI_INT, MPI_SUM, comm.comm_);
#endif
    }

    if (numThreadsOnThisNode > hwTop.maxThreads())
    {
        std::string mesg = "WARNING: ";
        if (GMX_LIB_MPI)
        {
            mesg += formatString("On rank %d: o", rank);
        }
        else
        {
            mesg += "O";
        }
        mesg += formatString("versubscribing the recommended max load of %d logical CPUs",
                             hwTop.maxThreads());
        if (GMX_LIB_MPI)
        {
            mesg += " per node";
        }
        mesg += formatString(" with %d ", numThreadsOnThisNode);
        if (numRanksOnThisNode == numThreadsOnThisNode)
        {
            if (GMX_THREAD_MPI)
            {
                mesg += "thread-MPI threads.";
            }
            else
            {
                mesg += "MPI processes.";
            }
        }
        else
        {
            mesg += "threads.";
        }
        mesg += "\n         This will cause considerable performance loss.";
        /* Note that only the main rank logs to stderr and only ranks
         * with an open log file write to log.
         * TODO: When we have a proper parallel logging framework,
         *       the framework should add the rank and node numbers.
         */
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted("%s", mesg.c_str());
    }
}

} // namespace gmx
