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

#include "gmx_omp_nthreads.h"

#include "config.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <filesystem>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/programcontext.h"

/** Structure with the number of threads for each OpenMP multi-threaded
 *  algorithmic module in mdrun. */
typedef struct
{
    int gnth;     /**< Global num. of threads per PP or PP+PME process/tMPI thread. */
    int gnth_pme; /**< Global num. of threads per PME only process/tMPI thread. */

    gmx::EnumerationArray<ModuleMultiThread, int> nth; /**< Number of threads for each module, indexed with module_nth_t */
} omp_module_nthreads_t;

/** Names of environment variables to set the per module number of threads.
 *
 *  Indexed with the values of module_nth_t.
 * */
static const char* enumValueToEnvVariableString(ModuleMultiThread enumValue)
{
    constexpr gmx::EnumerationArray<ModuleMultiThread, const char*> moduleMultiThreadEnvVariableNames = {
        "GMX_DEFAULT_NUM_THREADS should never be set",
        "GMX_DOMDEC_NUM_THREADS",
        "GMX_PAIRSEARCH_NUM_THREADS",
        "GMX_NONBONDED_NUM_THREADS",
        "GMX_LISTED_FORCES_NUM_THREADS",
        "GMX_PME_NUM_THREADS",
        "GMX_UPDATE_NUM_THREADS",
        "GMX_VSITE_NUM_THREADS",
        "GMX_LINCS_NUM_THREADS",
        "GMX_SETTLE_NUM_THREADS"
    };
    return moduleMultiThreadEnvVariableNames[enumValue];
}

/** Names of the modules. */
static const char* enumValueToString(ModuleMultiThread enumValue)
{
    constexpr gmx::EnumerationArray<ModuleMultiThread, const char*> moduleMultiThreadNames = {
        "default", "domain decomposition", "pair search", "non-bonded", "bonded", "PME",
        "update",  "virtual sites",        "LINCS",       "SETTLE"
    };
    return moduleMultiThreadNames[enumValue];
}

/** Number of threads for each algorithmic module.
 *
 *  File-scope global variable that gets set once in pick_module_nthreads()
 *  and queried via gmx_omp_nthreads_get().
 *
 *  All fields are initialized to 0 which should result in errors if
 *  the init call is omitted.
 * */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static omp_module_nthreads_t modth = { 0, 0, { 0, 0, 0, 0, 0, 0, 0, 0, 0 } };


/** Determine the number of threads for module \p mod.
 *
 *  \p m takes values form the module_nth_t enum and maps these to the
 *  corresponding value in modth_env_var.
 *
 *  Each number of threads per module takes the default value unless
 *  GMX_*_NUM_THERADS env var is set, case in which its value overrides
 *  the default.
 */
static void pick_module_nthreads(const gmx::MDLogger& mdlog, ModuleMultiThread m, gmx_bool bSepPME)
{
    char* env;
    int   nth;

    const bool bOMP = GMX_OPENMP;

    /* The default should never be set through a GMX_*_NUM_THREADS env var
     * as it's always equal with gnth. */
    if (m == ModuleMultiThread::Default)
    {
        return;
    }

    /* check the environment variable */
    if ((env = getenv(enumValueToEnvVariableString(m))) != nullptr)
    {
        sscanf(env, "%d", &nth);

        if (!bOMP)
        {
            gmx_warning("%s=%d is set, but %s is compiled without OpenMP!",
                        enumValueToEnvVariableString(m),
                        nth,
                        gmx::getProgramContext().displayName());
        }

        /* with the verlet codepath, when any GMX_*_NUM_THREADS env var is set,
         * OMP_NUM_THREADS also has to be set */
        if (getenv("OMP_NUM_THREADS") == nullptr)
        {
            gmx_warning(
                    "%s=%d is set, the default number of threads also "
                    "needs to be set with OMP_NUM_THREADS!",
                    enumValueToEnvVariableString(m),
                    nth);
        }

        /* only babble if we are really overriding with a different value */
        if ((bSepPME && m == ModuleMultiThread::Pme && nth != modth.gnth_pme) || (nth != modth.gnth))
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted("%s=%d set, overriding the default number of %s threads",
                                         enumValueToEnvVariableString(m),
                                         nth,
                                         enumValueToString(m));
        }
    }
    else
    {
        /* pick the global PME node nthreads if we are setting the number
         * of threads in separate PME nodes  */
        nth = (bSepPME && m == ModuleMultiThread::Pme) ? modth.gnth_pme : modth.gnth;
    }

    gmx_omp_nthreads_set(m, nth);
}

void gmx_omp_nthreads_read_env(const gmx::MDLogger& mdlog, int* nthreads_omp)
{
    char*    env;
    gmx_bool bCommandLineSetNthreadsOMP = *nthreads_omp > 0;
    char     buffer[STRLEN];

    GMX_RELEASE_ASSERT(nthreads_omp, "nthreads_omp must be a non-NULL pointer");

    if ((env = getenv("OMP_NUM_THREADS")) != nullptr)
    {
        int nt_omp;

        sscanf(env, "%d", &nt_omp);
        if (nt_omp <= 0)
        {
            gmx_fatal(FARGS, "OMP_NUM_THREADS is invalid: '%s'", env);
        }

        if (bCommandLineSetNthreadsOMP && nt_omp != *nthreads_omp)
        {
            gmx_fatal(FARGS,
                      "Environment variable OMP_NUM_THREADS (%d) and the number of threads "
                      "requested on the command line (%d) have different values. Either omit one, "
                      "or set them both to the same value.",
                      nt_omp,
                      *nthreads_omp);
        }

        /* Setting the number of OpenMP threads. */
        *nthreads_omp = nt_omp;

        /* Output the results */
        sprintf(buffer,
                "\nThe number of OpenMP threads was set by environment variable OMP_NUM_THREADS to "
                "%d%s\n\n",
                nt_omp,
                bCommandLineSetNthreadsOMP ? " (and the command-line setting agreed with that)" : "");

        /* This prints once per simulation for multi-simulations,
         * which might help diagnose issues with inhomogenous
         * cluster setups. */
        GMX_LOG(mdlog.info).appendTextFormatted("%s", buffer);
        if (debug)
        {
            /* This prints once per process for real MPI (i.e. once
             * per debug file), and once per simulation for thread MPI
             * (because of logic in the calling function). */
            fputs(buffer, debug);
        }
    }
}

/*! \brief Helper function for parsing various input about the number
    of OpenMP threads to use in various modules and deciding what to
    do about it. */
static void manage_number_of_openmp_threads(const gmx::MDLogger& mdlog,
                                            const t_commrec*     cr,
                                            bool                 bOMP,
                                            int                  maxThreads,
                                            int                  omp_nthreads_req,
                                            int                  omp_nthreads_pme_req,
                                            gmx_bool gmx_unused  bThisNodePMEOnly,
                                            int                  numRanksOnThisNode,
                                            gmx_bool             bSepPME)
{
    int   nth;
    char* env;
    bool  threadLimitApplied{ false };

#if GMX_THREAD_MPI
    /* modth is shared among tMPI threads, so for thread safety, the
     * detection is done on the main only. It is not thread-safe
     * with multiple simulations, but that's anyway not supported by
     * tMPI. */
    if (!SIMMAIN(cr))
    {
        return;
    }
#else
    GMX_UNUSED_VALUE(cr);
#endif

    /* With full OpenMP support (verlet scheme) set the number of threads
     * per process / default:
     * - 1 if not compiled with OpenMP or
     * - OMP_NUM_THREADS if the env. var is set, or
     * - omp_nthreads_req = #of threads requested by the user on the mdrun
     *   command line, otherwise
     * - take the max number of available threads and distribute them
     *   on the processes/tMPI threads.
     * ~ The GMX_*_NUM_THREADS env var overrides the number of threads of
     *   the respective module and it has to be used in conjunction with
     *   OMP_NUM_THREADS.
     *
     * The number of PME threads is equal to:
     * - 1 if not compiled with OpenMP or
     * - GMX_PME_NUM_THREADS if defined, otherwise
     * - OMP_NUM_THREADS if defined, otherwise
     * - 1
     */
    nth = 1;
    if ((env = getenv("OMP_NUM_THREADS")) != nullptr)
    {
        if (!bOMP && (std::strncmp(env, "1", 1) != 0))
        {
            gmx_warning("OMP_NUM_THREADS is set, but %s was compiled without OpenMP support!",
                        gmx::getProgramContext().displayName());
        }
        else
        {
            nth = gmx_omp_get_max_threads();
        }
    }
    else if (omp_nthreads_req > 0)
    {
        nth = omp_nthreads_req;
    }
    else if (bOMP)
    {
        /* max available threads per node */
        nth = maxThreads;

        /* divide the threads among the MPI ranks */
        if (nth >= numRanksOnThisNode)
        {
            nth /= numRanksOnThisNode;
        }
        else
        {
            nth = 1;
        }
    }

    if (nth > GMX_OPENMP_MAX_THREADS)
    {
        nth                = GMX_OPENMP_MAX_THREADS;
        threadLimitApplied = true;
    }

    /* now we have the global values, set them:
     * - 1 if not compiled with OpenMP
     * - nth for the verlet scheme when compiled with OpenMP
     */
    if (bOMP)
    {
        modth.gnth = nth;
    }
    else
    {
        modth.gnth = 1;
    }

    if (bSepPME)
    {
        if (omp_nthreads_pme_req > 0)
        {
            modth.gnth_pme = omp_nthreads_pme_req;
        }
        else
        {
            modth.gnth_pme = nth;
        }
    }
    else
    {
        modth.gnth_pme = 0;
    }

    if (modth.gnth_pme > GMX_OPENMP_MAX_THREADS)
    {
        modth.gnth_pme     = GMX_OPENMP_MAX_THREADS;
        threadLimitApplied = true;
    }

    if (threadLimitApplied)
    {
        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "Applying OpenMP thread count limit of %d (imposed by the "
                        "GMX_OPENMP_MAX_THREADS compile-time setting).",
                        GMX_OPENMP_MAX_THREADS);
    }

    /* now set the per-module values */
    modth.nth[ModuleMultiThread::Default] = modth.gnth;
    pick_module_nthreads(mdlog, ModuleMultiThread::Domdec, bSepPME);
    pick_module_nthreads(mdlog, ModuleMultiThread::Pairsearch, bSepPME);
    pick_module_nthreads(mdlog, ModuleMultiThread::Nonbonded, bSepPME);
    pick_module_nthreads(mdlog, ModuleMultiThread::Bonded, bSepPME);
    pick_module_nthreads(mdlog, ModuleMultiThread::Pme, bSepPME);
    pick_module_nthreads(mdlog, ModuleMultiThread::Update, bSepPME);
    pick_module_nthreads(mdlog, ModuleMultiThread::VirtualSite, bSepPME);
    pick_module_nthreads(mdlog, ModuleMultiThread::Lincs, bSepPME);
    pick_module_nthreads(mdlog, ModuleMultiThread::Settle, bSepPME);

    /* set the number of threads globally */
    if (bOMP)
    {
#if !GMX_THREAD_MPI
        if (bThisNodePMEOnly)
        {
            gmx_omp_set_num_threads(modth.gnth_pme);
        }
        else
#endif /* GMX_THREAD_MPI */
        {
            gmx_omp_set_num_threads(nth);
        }
    }
}

/*! \brief Report on the OpenMP settings that will be used */
static void reportOpenmpSettings(const gmx::MDLogger& mdlog, const t_commrec* cr, gmx_bool bOMP, gmx_bool bSepPME)
{
#if GMX_THREAD_MPI
    const char* mpi_str = "per tMPI thread";
#else
    const char* mpi_str = "per MPI process";
#endif
    int nth_min, nth_max, nth_pme_min, nth_pme_max;

    /* inform the user about the settings */
    if (!bOMP)
    {
        return;
    }

#if GMX_MPI
    if (cr->nnodes > 1)
    {
        /* Get the min and max thread counts over the MPI ranks */
        int buf_in[4], buf_out[4];

        buf_in[0] = -modth.gnth;
        buf_in[1] = modth.gnth;
        buf_in[2] = -modth.gnth_pme;
        buf_in[3] = modth.gnth_pme;

        MPI_Allreduce(buf_in, buf_out, 4, MPI_INT, MPI_MAX, cr->mpi_comm_mysim);

        nth_min     = -buf_out[0];
        nth_max     = buf_out[1];
        nth_pme_min = -buf_out[2];
        nth_pme_max = buf_out[3];
    }
    else
#endif
    {
        nth_min     = modth.gnth;
        nth_max     = modth.gnth;
        nth_pme_min = modth.gnth_pme;
        nth_pme_max = modth.gnth_pme;
    }


    if (nth_max == nth_min)
    {
        GMX_LOG(mdlog.warning)
                .appendTextFormatted("Using %d OpenMP thread%s %s",
                                     nth_min,
                                     nth_min > 1 ? "s" : "",
                                     cr->nnodes > 1 ? mpi_str : "");
    }
    else
    {
        GMX_LOG(mdlog.warning).appendTextFormatted("Using %d - %d OpenMP threads %s", nth_min, nth_max, mpi_str);
    }

    if (bSepPME && (nth_pme_min != nth_min || nth_pme_max != nth_max))
    {
        if (nth_pme_max == nth_pme_min)
        {
            GMX_LOG(mdlog.warning)
                    .appendTextFormatted("Using %d OpenMP thread%s %s for PME",
                                         nth_pme_min,
                                         nth_pme_min > 1 ? "s" : "",
                                         cr->nnodes > 1 ? mpi_str : "");
        }
        else
        {
            GMX_LOG(mdlog.warning)
                    .appendTextFormatted(
                            "Using %d - %d OpenMP threads %s for PME", nth_pme_min, nth_pme_max, mpi_str);
        }
    }
    GMX_LOG(mdlog.warning);
}

void gmx_omp_nthreads_init(const gmx::MDLogger& mdlog,
                           t_commrec*           cr,
                           int                  maxThreads,
                           int                  numRanksOnThisNode,
                           int                  omp_nthreads_req,
                           int                  omp_nthreads_pme_req,
                           gmx_bool             bThisNodePMEOnly)
{
    gmx_bool bSepPME;

    const bool bOMP = GMX_OPENMP;

    bSepPME = (thisRankHasDuty(cr, DUTY_PP) != thisRankHasDuty(cr, DUTY_PME));

    manage_number_of_openmp_threads(
            mdlog, cr, bOMP, maxThreads, omp_nthreads_req, omp_nthreads_pme_req, bThisNodePMEOnly, numRanksOnThisNode, bSepPME);
#if GMX_THREAD_MPI
    /* Non-main threads have to wait for the OpenMP management to be
     * done, so that code elsewhere that uses OpenMP can be certain
     * the setup is complete. */
    if (PAR(cr))
    {
        MPI_Barrier(cr->mpi_comm_mysim);
    }
#endif

    reportOpenmpSettings(mdlog, cr, bOMP, bSepPME);
}

int gmx_omp_nthreads_get(ModuleMultiThread mod)
{
    if (mod < ModuleMultiThread::Default || mod >= ModuleMultiThread::Count)
    {
        /* invalid module queried */
        return -1;
    }
    else
    {
        return modth.nth[mod];
    }
}

void gmx_omp_nthreads_set(ModuleMultiThread mod, int nthreads)
{
    /* Catch an attempt to set the number of threads on an invalid
     * OpenMP module. */
    GMX_RELEASE_ASSERT(mod >= ModuleMultiThread::Default && mod < ModuleMultiThread::Count,
                       "Trying to set nthreads on invalid OpenMP module");

    modth.nth[mod] = nthreads;
}
