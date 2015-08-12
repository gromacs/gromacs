/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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

#include "resource-division.h"

#include "config.h"

#include <stdlib.h>
#include <string.h>

#include <algorithm>

#include "gromacs/legacyheaders/gmx_detect_hardware.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/md_logging.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"


/* DISCLAIMER: All the atom count and thread numbers below are heuristic.
 * The real switching points will depend on the system simulation,
 * the algorithms used and the hardware it's running on, as well as if there
 * are other jobs running on the same machine. We try to take into account
 * factors that have a large influence, such as recent Intel CPUs being
 * much better at wide multi-threading. The remaining factors should
 * (hopefully) have a small influence, such that the performance just before
 * and after a switch point doesn't change too much.
 */

#ifdef GMX_OPENMP
static const bool bHasOmpSupport = true;
#else
static const bool bHasOmpSupport = false;
#endif

#ifdef GMX_THREAD_MPI
/* The minimum number of atoms per tMPI thread. With fewer atoms than this,
 * the number of threads will get lowered.
 */
static const int min_atoms_per_mpi_thread =  90;
static const int min_atoms_per_gpu        = 900;
#endif /* GMX_THREAD_MPI */

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
#if defined GMX_OPENMP && defined GMX_MPI
const int nthreads_omp_mpi_ok_max              =  8;
const int nthreads_omp_mpi_ok_min_cpu          =  1;
#endif
const int nthreads_omp_mpi_ok_min_gpu          =  2;
const int nthreads_omp_mpi_target_max          =  6;


/* Returns the maximum OpenMP thread count for which using a single MPI rank
 * should be faster than using multiple ranks with the same total thread count.
 */
static int nthreads_omp_faster(gmx_cpuid_t cpuid_info, gmx_bool bUseGPU)
{
    int nth;

    if (gmx_cpuid_vendor(cpuid_info) == GMX_CPUID_VENDOR_INTEL &&
        gmx_cpuid_feature(cpuid_info, GMX_CPUID_FEATURE_X86_AVX))
    {
        nth = nthreads_omp_faster_Intel_AVX;
    }
    else if (gmx_cpuid_is_intel_nehalem(cpuid_info))
    {
        nth = nthreads_omp_faster_Nehalem;
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

/* Returns that maximum OpenMP thread count that passes the efficiency check */
static int nthreads_omp_efficient_max(int gmx_unused nrank,
                                      gmx_cpuid_t    cpuid_info,
                                      gmx_bool       bUseGPU)
{
#if defined GMX_OPENMP && defined GMX_MPI
    if (nrank > 1)
    {
        return nthreads_omp_mpi_ok_max;
    }
    else
#endif
    {
        return nthreads_omp_faster(cpuid_info, bUseGPU);
    }
}

/* Return the number of thread-MPI ranks to use.
 * This is chosen such that we can always obey our own efficiency checks.
 */
static int get_tmpi_omp_thread_division(const gmx_hw_info_t *hwinfo,
                                        const gmx_hw_opt_t  *hw_opt,
                                        int                  nthreads_tot,
                                        int                  ngpu)
{
    int nrank;

    GMX_RELEASE_ASSERT(nthreads_tot > 0, "There must be at least one thread per rank");

    /* There are no separate PME nodes here, as we ensured in
     * check_and_update_hw_opt that nthreads_tmpi>0 with PME nodes
     * and a conditional ensures we would not have ended up here.
     * Note that separate PME nodes might be switched on later.
     */
    if (ngpu > 0)
    {
        nrank = ngpu;
        if (nthreads_tot < nrank)
        {
            /* #thread < #gpu is very unlikely, but if so: waste gpu(s) */
            nrank = nthreads_tot;
        }
        else if (gmx_gpu_sharing_supported() &&
                 (nthreads_tot > nthreads_omp_faster(hwinfo->cpuid_info,
                                                     ngpu > 0) ||
                  (ngpu > 1 && nthreads_tot/ngpu > nthreads_omp_mpi_target_max)))
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
    else if (hw_opt->nthreads_omp > 0)
    {
        /* Here we could oversubscribe, when we do, we issue a warning later */
        nrank = std::max(1, nthreads_tot/hw_opt->nthreads_omp);
    }
    else
    {
        if (nthreads_tot <= nthreads_omp_faster(hwinfo->cpuid_info, ngpu > 0))
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


static int getMaxGpuUsable(FILE *fplog, const t_commrec *cr, const gmx_hw_info_t *hwinfo, int cutoff_scheme)
{
    /* This code relies on the fact that GPU are not detected when GPU
     * acceleration was disabled at run time by the user.
     */
    if (cutoff_scheme == ecutsVERLET &&
        hwinfo->gpu_info.n_dev_compatible > 0)
    {
        if (gmx_multiple_gpu_per_node_supported())
        {
            return hwinfo->gpu_info.n_dev_compatible;
        }
        else
        {
            if (hwinfo->gpu_info.n_dev_compatible > 1)
            {
                md_print_warn(cr, fplog, "More than one compatible GPU is available, but GROMACS can only use one of them. Using a single thread-MPI rank.\n");
            }
            return 1;
        }
    }
    else
    {
        return 0;
    }
}


#ifdef GMX_THREAD_MPI
/* Get the number of MPI ranks to use for thread-MPI based on how many
 * were requested, which algorithms we're using,
 * and how many particles there are.
 * At the point we have already called check_and_update_hw_opt.
 * Thus all options should be internally consistent and consistent
 * with the hardware, except that ntmpi could be larger than #GPU.
 */
int get_nthreads_mpi(const gmx_hw_info_t *hwinfo,
                     gmx_hw_opt_t        *hw_opt,
                     const t_inputrec    *inputrec,
                     const gmx_mtop_t    *mtop,
                     const t_commrec     *cr,
                     FILE                *fplog)
{
    int      nthreads_hw, nthreads_tot_max, nrank, ngpu;
    int      min_atoms_per_mpi_rank;

    /* Check if an algorithm does not support parallel simulation.  */
    if (inputrec->eI == eiLBFGS ||
        inputrec->coulombtype == eelEWALD)
    {
        md_print_warn(cr, fplog, "The integration or electrostatics algorithm doesn't support parallel runs. Using a single thread-MPI rank.\n");
        if (hw_opt->nthreads_tmpi > 1)
        {
            gmx_fatal(FARGS, "You asked for more than 1 thread-MPI rank, but an algorithm doesn't support that");
        }

        return 1;
    }

    if (hw_opt->nthreads_tmpi > 0)
    {
        /* Trivial, return right away */
        return hw_opt->nthreads_tmpi;
    }

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

    ngpu = getMaxGpuUsable(fplog, cr, hwinfo, inputrec->cutoff_scheme);

    if (inputrec->cutoff_scheme == ecutsGROUP)
    {
        /* We checked this before, but it doesn't hurt to do it once more */
        GMX_RELEASE_ASSERT(hw_opt->nthreads_omp == 1, "The group scheme only supports one OpenMP thread per rank");
    }

    nrank =
        get_tmpi_omp_thread_division(hwinfo, hw_opt, nthreads_tot_max, ngpu);

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
        if (gmx_cpuid_x86_smt(hwinfo->cpuid_info) == GMX_CPUID_X86_SMT_ENABLED &&
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

            nt_omp_max = nthreads_omp_efficient_max(nrank, hwinfo->cpuid_info, ngpu >= 1);

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
                                        const gmx_hw_opt_t  *hw_opt,
                                        gmx_bool             bNtOmpOptionSet,
                                        t_commrec           *cr,
                                        FILE                *fplog)
{
#if defined GMX_OPENMP && defined GMX_MPI
    int         nth_omp_min, nth_omp_max, ngpu;
    char        buf[1000];
#ifdef GMX_THREAD_MPI
    const char *mpi_option = " (option -ntmpi)";
#else
    const char *mpi_option = "";
#endif

    /* This function should be called after thread-MPI (when configured) and
     * OpenMP have been initialized. Check that here.
     */
#ifdef GMX_THREAD_MPI
    GMX_RELEASE_ASSERT(nthreads_omp_faster_default >= nthreads_omp_mpi_ok_max, "Inconsistent OpenMP thread count default values");
    GMX_RELEASE_ASSERT(hw_opt->nthreads_tmpi >= 1, "Must have at least one thread-MPI rank");
#endif
    GMX_RELEASE_ASSERT(gmx_omp_nthreads_get(emntDefault) >= 1, "Must have at least one OpenMP thread");

    nth_omp_min = gmx_omp_nthreads_get(emntDefault);
    nth_omp_max = gmx_omp_nthreads_get(emntDefault);
    ngpu        = hw_opt->gpu_opt.n_dev_use;

    /* Thread-MPI seems to have a bug with reduce on 1 node, so use a cond. */
    if (cr->nnodes + cr->npmenodes > 1)
    {
        int count[3], count_max[3];

        count[0] = -nth_omp_min;
        count[1] =  nth_omp_max;
        count[2] =  ngpu;

        MPI_Allreduce(count, count_max, 3, MPI_INT, MPI_MAX, cr->mpi_comm_mysim);

        /* In case of an inhomogeneous run setup we use the maximum counts */
        nth_omp_min = -count_max[0];
        nth_omp_max =  count_max[1];
        ngpu        =  count_max[2];
    }

    int nthreads_omp_mpi_ok_min;

    if (ngpu == 0)
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
            (!(ngpu > 0 && !gmx_gpu_sharing_supported()) &&
             nth_omp_max > nthreads_omp_mpi_ok_max))
        {
            /* Note that we print target_max here, not ok_max */
            sprintf(buf, "Your choice of number of MPI ranks and amount of resources results in using %d OpenMP threads per rank, which is most likely inefficient. The optimum is usually between %d and %d threads per rank.",
                    nth_omp_max,
                    nthreads_omp_mpi_ok_min,
                    nthreads_omp_mpi_target_max);

            if (bNtOmpOptionSet)
            {
                md_print_warn(cr, fplog, "NOTE: %s\n", buf);
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
        /* No domain decomposition (or only one domain) */
        if (!(ngpu > 0 && !gmx_gpu_sharing_supported()) &&
            nth_omp_max > nthreads_omp_faster(hwinfo->cpuid_info, ngpu > 0))
        {
            /* To arrive here, the user/system set #ranks and/or #OMPthreads */
            gmx_bool bEnvSet;
            char     buf2[256];

            bEnvSet = (getenv("OMP_NUM_THREADS") != NULL);

            if (bNtOmpOptionSet || bEnvSet)
            {
                sprintf(buf2, "You requested %d OpenMP threads", nth_omp_max);
            }
            else
            {
                sprintf(buf2, "Your choice of %d MPI rank%s and the use of %d total threads %sleads to the use of %d OpenMP threads",
                        cr->nnodes + cr->npmenodes,
                        cr->nnodes + cr->npmenodes == 1 ? "" : "s",
                        hw_opt->nthreads_tot > 0 ? hw_opt->nthreads_tot : hwinfo->nthreads_hw_avail,
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
                md_print_warn(cr, fplog, "NOTE: %s\n", buf);
            }
            else
            {
                gmx_fatal(FARGS, "%s If you want to run with this many OpenMP threads, specify the -ntomp option. But we suggest to increase the number of MPI ranks%s.", buf, mpi_option);
            }
        }
    }
#else /* GMX_OPENMP && GMX_MPI */
      /* No OpenMP and/or MPI: it doesn't make much sense to check */
    GMX_UNUSED_VALUE(hw_opt);
    GMX_UNUSED_VALUE(bNtOmpOptionSet);
    /* Check if we have more than 1 physical core, if detected,
     * or more than 1 hardware thread if physical cores were not detected.
     */
#if !(defined GMX_OPENMP) && !(defined GMX_MPI)
    if ((hwinfo->ncore > 1) ||
        (hwinfo->ncore == 0 && hwinfo->nthreads_hw_avail > 1))
    {
        md_print_warn(cr, fplog, "NOTE: GROMACS was compiled without OpenMP and (thread-)MPI support, can only use a single CPU core\n");
    }
#else
    GMX_UNUSED_VALUE(hwinfo);
    GMX_UNUSED_VALUE(cr);
    GMX_UNUSED_VALUE(fplog);
#endif

#endif /* GMX_OPENMP && GMX_MPI */
}


static void print_hw_opt(FILE *fp, const gmx_hw_opt_t *hw_opt)
{
    fprintf(fp, "hw_opt: nt %d ntmpi %d ntomp %d ntomp_pme %d gpu_id '%s'\n",
            hw_opt->nthreads_tot,
            hw_opt->nthreads_tmpi,
            hw_opt->nthreads_omp,
            hw_opt->nthreads_omp_pme,
            hw_opt->gpu_opt.gpu_id != NULL ? hw_opt->gpu_opt.gpu_id : "");
}

/* Checks we can do when we don't (yet) know the cut-off scheme */
void check_and_update_hw_opt_1(gmx_hw_opt_t    *hw_opt,
                               const t_commrec *cr)
{
    /* Currently hw_opt only contains default settings or settings supplied
     * by the user on the command line.
     * Check for OpenMP settings stored in environment variables, which can
     * potentially be different on different MPI ranks.
     */
    gmx_omp_nthreads_read_env(&hw_opt->nthreads_omp, SIMMASTER(cr));

    /* Check restrictions on the user supplied options before modifying them.
     * TODO: Put the user values in a const struct and preserve them.
     */
#ifndef GMX_THREAD_MPI
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
            cr->npmenodes <= 0)
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

    if (hw_opt->nthreads_tot == 1)
    {
        hw_opt->nthreads_tmpi = 1;

        if (hw_opt->nthreads_omp > 1)
        {
            gmx_fatal(FARGS, "You requested %d OpenMP threads with %d total threads",
                      hw_opt->nthreads_tmpi, hw_opt->nthreads_tot);
        }
        hw_opt->nthreads_omp     = 1;
        hw_opt->nthreads_omp_pme = 1;
    }

    /* Parse GPU IDs, if provided.
     * We check consistency with the tMPI thread count later.
     */
    gmx_parse_gpu_ids(&hw_opt->gpu_opt);

#ifdef GMX_THREAD_MPI
    if (hw_opt->gpu_opt.n_dev_use > 0 && hw_opt->nthreads_tmpi == 0)
    {
        /* Set the number of MPI threads equal to the number of GPUs */
        hw_opt->nthreads_tmpi = hw_opt->gpu_opt.n_dev_use;

        if (hw_opt->nthreads_tot > 0 &&
            hw_opt->nthreads_tmpi > hw_opt->nthreads_tot)
        {
            /* We have more GPUs than total threads requested.
             * We choose to (later) generate a mismatch error,
             * instead of launching more threads than requested.
             */
            hw_opt->nthreads_tmpi = hw_opt->nthreads_tot;
        }
    }
#endif

    if (debug)
    {
        print_hw_opt(debug, hw_opt);
    }

    /* Asserting this simplifies the hardware resource division later
     * on. */
    GMX_RELEASE_ASSERT(!(hw_opt->nthreads_omp_pme >= 1 && hw_opt->nthreads_omp <= 0),
                       "PME thread count should only be set when the normal thread count is also set");
}

/* Checks we can do when we know the cut-off scheme */
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

/* Checks we can do when we know the thread-MPI rank count */
void check_and_update_hw_opt_3(gmx_hw_opt_t gmx_unused *hw_opt)
{
#ifdef GMX_THREAD_MPI
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
