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

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>

#include "gromacs/legacyheaders/gmx_detect_hardware.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/md_logging.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/utility/fatalerror.h"


#ifdef GMX_THREAD_MPI
/* The minimum number of atoms per tMPI thread. With fewer atoms than this,
 * the number of threads will get lowered.
 */
static const int min_atoms_per_mpi_thread =  90;
static const int min_atoms_per_gpu        = 900;

static int get_tmpi_omp_thread_division(const gmx_hw_info_t *hwinfo,
                                        const gmx_hw_opt_t  *hw_opt,
                                        int                  nthreads_tot,
                                        int                  ngpu)
{
    int nthreads_tmpi;

    /* There are no separate PME nodes here, as we ensured in
     * check_and_update_hw_opt that nthreads_tmpi>0 with PME nodes
     * and a conditional ensures we would not have ended up here.
     * Note that separate PME nodes might be switched on later.
     */
    if (ngpu > 0)
    {
        nthreads_tmpi = ngpu;

        /* When the user sets nthreads_omp, we can end up oversubscribing CPU cores
         * if we simply start as many ranks as GPUs. To avoid this, we start as few
         * tMPI ranks as necessary to avoid oversubscription and instead leave GPUs idle.
         * If the user does not set the number of OpenMP threads, nthreads_omp==0 and
         * this code has no effect.
         */
        while (nthreads_tmpi*hw_opt->nthreads_omp > hwinfo->nthreads_hw_avail && nthreads_tmpi > 1)
        {
            nthreads_tmpi--;
        }


        if (nthreads_tot > 0 && nthreads_tot < nthreads_tmpi)
        {
            nthreads_tmpi = nthreads_tot;
        }
    }
    else if (hw_opt->nthreads_omp > 0)
    {
        /* Here we could oversubscribe, when we do, we issue a warning later */
        nthreads_tmpi = std::max(1, nthreads_tot/hw_opt->nthreads_omp);
    }
    else
    {
        /* TODO choose nthreads_omp based on hardware topology
           when we have a hardware topology detection library */
        /* In general, when running up to 4 threads, OpenMP should be faster.
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
        const int nthreads_omp_always_faster             =  4;
        const int nthreads_omp_always_faster_Nehalem     = 12;
        const int nthreads_omp_always_faster_Intel_AVX   = 16;
        gmx_bool  bIntelAVX;

        bIntelAVX =
            (gmx_cpuid_vendor(hwinfo->cpuid_info) == GMX_CPUID_VENDOR_INTEL &&
             gmx_cpuid_feature(hwinfo->cpuid_info, GMX_CPUID_FEATURE_X86_AVX));

        if (nthreads_tot <= nthreads_omp_always_faster ||
            ((gmx_cpuid_is_intel_nehalem(hwinfo->cpuid_info) && nthreads_tot <= nthreads_omp_always_faster_Nehalem) ||
             (bIntelAVX && nthreads_tot <= nthreads_omp_always_faster_Intel_AVX)))
        {
            /* Use pure OpenMP parallelization */
            nthreads_tmpi = 1;
        }
        else
        {
            /* Don't use OpenMP parallelization */
            nthreads_tmpi = nthreads_tot;
        }
    }

    return nthreads_tmpi;
}


/* Get the number of threads to use for thread-MPI based on how many
 * were requested, which algorithms we're using,
 * and how many particles there are.
 * At the point we have already called check_and_update_hw_opt.
 * Thus all options should be internally consistent and consistent
 * with the hardware, except that ntmpi could be larger than #GPU.
 */
int get_nthreads_mpi(const gmx_hw_info_t *hwinfo,
                     const gmx_hw_opt_t  *hw_opt,
                     const t_inputrec    *inputrec,
                     const gmx_mtop_t    *mtop,
                     const t_commrec     *cr,
                     FILE                *fplog)
{
    int      nthreads_hw, nthreads_tot_max, nthreads_tmpi, nthreads_new, ngpu;
    int      min_atoms_per_mpi_rank;
    gmx_bool bCanUseGPU;

    if (hw_opt->nthreads_tmpi > 0)
    {
        /* Trivial, return right away */
        return hw_opt->nthreads_tmpi;
    }

    nthreads_hw = hwinfo->nthreads_hw_avail;

    /* How many total (#tMPI*#OpenMP) threads can we start? */
    if (hw_opt->nthreads_tot > 0)
    {
        nthreads_tot_max = hw_opt->nthreads_tot;
    }
    else
    {
        nthreads_tot_max = nthreads_hw;
    }

    bCanUseGPU = (inputrec->cutoff_scheme == ecutsVERLET &&
                  hwinfo->gpu_info.n_dev_compatible > 0);
    if (bCanUseGPU)
    {
        ngpu = hwinfo->gpu_info.n_dev_compatible;
    }
    else
    {
        ngpu = 0;
    }

    if (inputrec->cutoff_scheme == ecutsGROUP)
    {
        /* We checked this before, but it doesn't hurt to do it once more */
        assert(hw_opt->nthreads_omp == 1);
    }

    nthreads_tmpi =
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
        if (bCanUseGPU)
        {
            min_atoms_per_mpi_rank = min_atoms_per_gpu;
        }
        else
        {
            min_atoms_per_mpi_rank = min_atoms_per_mpi_thread;
        }
    }

    /* Check if an algorithm does not support parallel simulation.  */
    if (nthreads_tmpi != 1 &&
        ( inputrec->eI == eiLBFGS ||
          inputrec->coulombtype == eelEWALD ) )
    {
        nthreads_tmpi = 1;

        md_print_warn(cr, fplog, "The integration or electrostatics algorithm doesn't support parallel runs. Using a single thread-MPI thread.\n");
        if (hw_opt->nthreads_tmpi > nthreads_tmpi)
        {
            gmx_fatal(FARGS, "You asked for more than 1 thread-MPI thread, but an algorithm doesn't support that");
        }
    }
    else if (mtop->natoms/nthreads_tmpi < min_atoms_per_mpi_rank)
    {
        /* the thread number was chosen automatically, but there are too many
           threads (too few atoms per thread) */
        nthreads_new = std::max(1, mtop->natoms/min_atoms_per_mpi_rank);

        /* Avoid partial use of Hyper-Threading */
        if (gmx_cpuid_x86_smt(hwinfo->cpuid_info) == GMX_CPUID_X86_SMT_ENABLED &&
            nthreads_new > nthreads_hw/2 && nthreads_new < nthreads_hw)
        {
            nthreads_new = nthreads_hw/2;
        }

        /* Avoid large prime numbers in the thread count */
        if (nthreads_new >= 6)
        {
            /* Use only 6,8,10 with additional factors of 2 */
            int fac;

            fac = 2;
            while (3*fac*2 <= nthreads_new)
            {
                fac *= 2;
            }

            nthreads_new = (nthreads_new/fac)*fac;
        }
        else
        {
            /* Avoid 5 */
            if (nthreads_new == 5)
            {
                nthreads_new = 4;
            }
        }

        nthreads_tmpi = nthreads_new;

        fprintf(stderr, "\n");
        fprintf(stderr, "NOTE: Parallelization is limited by the small number of atoms,\n");
        fprintf(stderr, "      only starting %d thread-MPI threads.\n", nthreads_tmpi);
        fprintf(stderr, "      You can use the -nt and/or -ntmpi option to optimize the number of threads.\n\n");
    }

    return nthreads_tmpi;
}
#endif /* GMX_THREAD_MPI */


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
void check_and_update_hw_opt_1(gmx_hw_opt_t *hw_opt,
                               gmx_bool      bIsSimMaster)
{
    gmx_omp_nthreads_read_env(&hw_opt->nthreads_omp, bIsSimMaster);

#ifndef GMX_THREAD_MPI
    if (hw_opt->nthreads_tot > 0)
    {
        gmx_fatal(FARGS, "Setting the total number of threads is only supported with thread-MPI and GROMACS was compiled without thread-MPI");
    }
    if (hw_opt->nthreads_tmpi > 0)
    {
        gmx_fatal(FARGS, "Setting the number of thread-MPI threads is only supported with thread-MPI and GROMACS was compiled without thread-MPI");
    }
#endif

#ifndef GMX_OPENMP
    if (hw_opt->nthreads_omp > 1)
    {
        gmx_fatal(FARGS, "More than 1 OpenMP thread requested, but GROMACS was compiled without OpenMP support");
    }
    hw_opt->nthreads_omp = 1;
#endif

    if (hw_opt->nthreads_tot > 0 && hw_opt->nthreads_omp_pme <= 0)
    {
        /* We have the same number of OpenMP threads for PP and PME processes,
         * thus we can perform several consistency checks.
         */
        if (hw_opt->nthreads_tmpi > 0 &&
            hw_opt->nthreads_omp > 0 &&
            hw_opt->nthreads_tot != hw_opt->nthreads_tmpi*hw_opt->nthreads_omp)
        {
            gmx_fatal(FARGS, "The total number of threads requested (%d) does not match the thread-MPI threads (%d) times the OpenMP threads (%d) requested",
                      hw_opt->nthreads_tot, hw_opt->nthreads_tmpi, hw_opt->nthreads_omp);
        }

        if (hw_opt->nthreads_tmpi > 0 &&
            hw_opt->nthreads_tot % hw_opt->nthreads_tmpi != 0)
        {
            gmx_fatal(FARGS, "The total number of threads requested (%d) is not divisible by the number of thread-MPI threads requested (%d)",
                      hw_opt->nthreads_tot, hw_opt->nthreads_tmpi);
        }

        if (hw_opt->nthreads_omp > 0 &&
            hw_opt->nthreads_tot % hw_opt->nthreads_omp != 0)
        {
            gmx_fatal(FARGS, "The total number of threads requested (%d) is not divisible by the number of OpenMP threads requested (%d)",
                      hw_opt->nthreads_tot, hw_opt->nthreads_omp);
        }

        if (hw_opt->nthreads_tmpi > 0 &&
            hw_opt->nthreads_omp <= 0)
        {
            hw_opt->nthreads_omp = hw_opt->nthreads_tot/hw_opt->nthreads_tmpi;
        }
    }

#ifndef GMX_OPENMP
    if (hw_opt->nthreads_omp > 1)
    {
        gmx_fatal(FARGS, "OpenMP threads are requested, but GROMACS was compiled without OpenMP support");
    }
#endif

    if (hw_opt->nthreads_omp_pme > 0 && hw_opt->nthreads_omp <= 0)
    {
        gmx_fatal(FARGS, "You need to specify -ntomp in addition to -ntomp_pme");
    }

    if (hw_opt->nthreads_tot == 1)
    {
        hw_opt->nthreads_tmpi = 1;

        if (hw_opt->nthreads_omp > 1)
        {
            gmx_fatal(FARGS, "You requested %d OpenMP threads with %d total threads",
                      hw_opt->nthreads_tmpi, hw_opt->nthreads_tot);
        }
        hw_opt->nthreads_omp = 1;
    }

    if (hw_opt->nthreads_omp_pme <= 0 && hw_opt->nthreads_omp > 0)
    {
        hw_opt->nthreads_omp_pme = hw_opt->nthreads_omp;
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

    if (hw_opt->nthreads_omp_pme <= 0 && hw_opt->nthreads_omp > 0)
    {
        hw_opt->nthreads_omp_pme = hw_opt->nthreads_omp;
    }

    if (debug)
    {
        print_hw_opt(debug, hw_opt);
    }
}
