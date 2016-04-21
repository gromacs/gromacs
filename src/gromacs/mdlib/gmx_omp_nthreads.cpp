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

#include "gmx_omp_nthreads.h"

#include "config.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/sts/sts.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/programcontext.h"

/** Structure with the number of threads for each OpenMP multi-threaded
 *  algorithmic module in mdrun. */
typedef struct
{
    int      gnth;          /**< Global num. of threads per PP or PP+PME process/tMPI thread. */
    int      gnth_pme;      /**< Global num. of threads per PME only process/tMPI thread. */

    int      nth[emntNR];   /**< Number of threads for each module, indexed with module_nth_t */
    gmx_bool initialized;   /**< TRUE if the module as been initialized. */
} omp_module_nthreads_t;

/** Names of environment variables to set the per module number of threads.
 *
 *  Indexed with the values of module_nth_t.
 * */
static const char *modth_env_var[emntNR] =
{
    "GMX_DEFAULT_NUM_THREADS should never be set",
    "GMX_DOMDEC_NUM_THREADS", "GMX_PAIRSEARCH_NUM_THREADS",
    "GMX_NONBONDED_NUM_THREADS", "GMX_LISTED_FORCES_NUM_THREADS",
    "GMX_PME_NUM_THREADS", "GMX_UPDATE_NUM_THREADS",
    "GMX_VSITE_NUM_THREADS",
    "GMX_LINCS_NUM_THREADS", "GMX_SETTLE_NUM_THREADS"
};

/** Names of the modules. */
static const char *mod_name[emntNR] =
{
    "default", "domain decomposition", "pair search", "non-bonded",
    "bonded", "PME", "update", "LINCS", "SETTLE"
};

/** Number of threads for each algorithmic module.
 *
 *  File-scope global variable that gets set once in pick_module_nthreads()
 *  and queried via gmx_omp_nthreads_get().
 *
 *  All fields are initialized to 0 which should result in errors if
 *  the init call is omitted.
 * */
static omp_module_nthreads_t modth = { 0, 0, {0, 0, 0, 0, 0, 0, 0, 0, 0}, FALSE};


/** Determine the number of threads for module \p mod.
 *
 *  \p m takes values form the module_nth_t enum and maps these to the
 *  corresponding value in modth_env_var.
 *
 *  Each number of threads per module takes the default value unless
 *  GMX_*_NUM_THERADS env var is set, case in which its value overrides
 *  the deafult.
 *
 *  The "group" scheme supports OpenMP only in PME and in thise case all but
 *  the PME nthread values default to 1.
 */
static void pick_module_nthreads(const gmx::MDLogger &mdlog, int m,
                                 gmx_bool bFullOmpSupport,
                                 gmx_bool bSepPME)
{
    char      *env;
    int        nth;

    const bool bOMP = GMX_STS;

    /* The default should never be set through a GMX_*_NUM_THREADS env var
     * as it's always equal with gnth. */
    if (m == emntDefault)
    {
        return;
    }

    /* check the environment variable */
    if ((env = getenv(modth_env_var[m])) != NULL)
    {
        sscanf(env, "%d", &nth);

        if (!bOMP)
        {
            gmx_warning("%s=%d is set, but %s is compiled without OpenMP!",
                        modth_env_var[m], nth,
                        gmx::getProgramContext().displayName());
        }

        /* with the verlet codepath, when any GMX_*_NUM_THREADS env var is set,
         * GMX_STS_THREADS also has to be set */
        if (bFullOmpSupport && getenv("GMX_STS_THREADS") == NULL)
        {
            gmx_warning("%s=%d is set, the default number of threads also "
                        "needs to be set with GMX_STS_THREADS!",
                        modth_env_var[m], nth);
        }

        /* with the group scheme warn if any env var except PME is set */
        if (!bFullOmpSupport)
        {
            if (m != emntPME)
            {
                gmx_warning("%s=%d is set, but OpenMP multithreading is not "
                            "supported in %s!",
                            modth_env_var[m], nth, mod_name[m]);
                nth = 1;
            }
        }

        /* only babble if we are really overriding with a different value */
        if ((bSepPME && m == emntPME && nth != modth.gnth_pme) || (nth != modth.gnth))
        {
            GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                    "%s=%d set, overriding the default number of %s threads",
                    modth_env_var[m], nth, mod_name[m]);
        }
    }
    else
    {
        /* pick the global PME node nthreads if we are setting the number
         * of threads in separate PME nodes  */
        nth = (bSepPME && m == emntPME) ? modth.gnth_pme : modth.gnth;
    }

    gmx_omp_nthreads_set(m, nth);
}

void gmx_omp_nthreads_read_env(int     *nthreads_omp,
                               gmx_bool bIsSimMaster)
{
    char    *env;
    gmx_bool bCommandLineSetNthreadsOMP = *nthreads_omp > 0;
    char     buffer[STRLEN];

    GMX_RELEASE_ASSERT(nthreads_omp, "nthreads_omp must be a non-NULL pointer");

    if ((env = getenv("GMX_STS_THREADS")) != NULL)
    {
        int nt_omp;

        sscanf(env, "%d", &nt_omp);
        if (nt_omp <= 0)
        {
            gmx_fatal(FARGS, "GMX_STS_THREADS is invalid: '%s'", env);
        }

        if (bCommandLineSetNthreadsOMP && nt_omp != *nthreads_omp)
        {
            gmx_fatal(FARGS, "Environment variable GMX_STS_THREADS (%d) and the number of threads requested on the command line (%d) have different values. Either omit one, or set them both to the same value.", nt_omp, *nthreads_omp);
        }

        /* Setting the number of OpenMP threads. */
        *nthreads_omp = nt_omp;

        /* Output the results */
        sprintf(buffer,
                "The number of OpenMP threads was set by environment variable GMX_STS_THREADS to %d%s\n",
                nt_omp,
                bCommandLineSetNthreadsOMP ? " (and the command-line setting agreed with that)" : "");
        if (bIsSimMaster)
        {
            /* This prints once per simulation for multi-simulations,
             * which might help diagnose issues with inhomogenous
             * cluster setups. */
            fputs(buffer, stderr);
        }
        if (debug)
        {
            /* This prints once per process for real MPI (i.e. once
             * per debug file), and once per simulation for thread MPI
             * (because of logic in the calling function). */
            fputs(buffer, debug);
        }
    }
    else if (!bCommandLineSetNthreadsOMP) {
        gmx_fatal(FARGS, "Either GMX_STS_THREADS must be set or number of threads given on the command line with -ntomp");
    }
}

/*! \brief Helper function for parsing various input about the number
    of OpenMP threads to use in various modules and deciding what to
    do about it. */
static void manage_number_of_openmp_threads(const gmx::MDLogger &mdlog,
                                            const t_commrec     *cr,
                                            bool                 bOMP,
                                            int      gmx_unused  nthreads_hw_avail,
                                            int                  omp_nthreads_req,
                                            int                  omp_nthreads_pme_req,
                                            gmx_bool gmx_unused  bThisNodePMEOnly,
                                            gmx_bool             bFullOmpSupport,
                                            int      gmx_unused  nppn,
                                            gmx_bool             bSepPME)
{
#if GMX_THREAD_MPI
    /* modth is shared among tMPI threads, so for thread safety, the
     * detection is done on the master only. It is not thread-safe
     * with multiple simulations, but that's anyway not supported by
     * tMPI. */
    if (!SIMMASTER(cr))
    {
        return;
    }
#else
    GMX_UNUSED_VALUE(cr);
#endif

    if (modth.initialized)
    {
        /* Just return if the initialization has already been
           done. This could only happen if gmx_omp_nthreads_init() has
           already been called. */
        return;
    }

    /* With full OpenMP support (verlet scheme) set the number of threads
     * per process / default:
     * - 1 if not compiled with OpenMP or
     * - GMX_STS_THREADS if the env. var is set, or
     * - omp_nthreads_req = #of threads requested by the user on the mdrun
     *   command line, otherwise
     * - take the max number of available threads and distribute them
     *   on the processes/tMPI threads.
     * ~ The GMX_*_NUM_THREADS env var overrides the number of threads of
     *   the respective module and it has to be used in conjunction with
     *   GMX_STS_THREADS.
     *
     * With the group scheme OpenMP multithreading is only supported in PME,
     * for all other modules nthreads is set to 1.
     * The number of PME threads is equal to:
     * - 1 if not compiled with OpenMP or
     * - GMX_PME_NUM_THREADS if defined, otherwise
     * - GMX_STS_THREADS if defined, otherwise
     * - 1
     */
    // Simplify greatly for STS
    int nth = omp_nthreads_req;
    gmx_omp_nthreads_read_env(&nth, true);
    STS::startup(nth);

    /* now we have the global values, set them:
     * - 1 if not compiled with OpenMP and for the group scheme
     * - nth for the verlet scheme when compiled with OpenMP
     */
    if (bFullOmpSupport && bOMP)
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

    /* now set the per-module values */
    modth.nth[emntDefault] = modth.gnth;
    pick_module_nthreads(mdlog, emntDomdec, bFullOmpSupport, bSepPME);
    pick_module_nthreads(mdlog, emntPairsearch, bFullOmpSupport, bSepPME);
    pick_module_nthreads(mdlog, emntNonbonded, bFullOmpSupport, bSepPME);
    pick_module_nthreads(mdlog, emntBonded, bFullOmpSupport, bSepPME);
    pick_module_nthreads(mdlog, emntPME, bFullOmpSupport, bSepPME);
    pick_module_nthreads(mdlog, emntUpdate, bFullOmpSupport, bSepPME);
    pick_module_nthreads(mdlog, emntVSITE, bFullOmpSupport, bSepPME);
    pick_module_nthreads(mdlog, emntLINCS, bFullOmpSupport, bSepPME);
    pick_module_nthreads(mdlog, emntSETTLE, bFullOmpSupport, bSepPME);

    /* set the number of threads globally */
    if (bOMP)
    {
#if !GMX_THREAD_MPI
        if (bThisNodePMEOnly)
        {
            gmx_omp_set_num_threads(modth.gnth_pme);
        }
        else
#endif      /* GMX_THREAD_MPI */
        {
            if (bFullOmpSupport)
            {
                gmx_omp_set_num_threads(nth);
            }
            else
            {
                gmx_omp_set_num_threads(1);
            }
        }
    }

    modth.initialized = TRUE;
}

/*! \brief Report on the OpenMP settings that will be used */
static void
reportOpenmpSettings(const gmx::MDLogger &mdlog,
                     const t_commrec     *cr,
                     gmx_bool             bOMP,
                     gmx_bool             bFullOmpSupport,
                     gmx_bool             bSepPME)
{
#if GMX_THREAD_MPI
    const char *mpi_str = "per tMPI thread";
#else
    const char *mpi_str = "per MPI process";
#endif
    int         nth_min, nth_max, nth_pme_min, nth_pme_max;

    /* inform the user about the settings */
    if (!bOMP)
    {
        return;
    }

#if GMX_MPI
    if (cr->nnodes + cr->npmenodes > 1)
    {
        /* Get the min and max thread counts over the MPI ranks */
        int buf_in[4], buf_out[4];

        buf_in[0] = -modth.gnth;
        buf_in[1] =  modth.gnth;
        buf_in[2] = -modth.gnth_pme;
        buf_in[3] =  modth.gnth_pme;

        MPI_Allreduce(buf_in, buf_out, 4, MPI_INT, MPI_MAX, cr->mpi_comm_mysim);

        nth_min     = -buf_out[0];
        nth_max     =  buf_out[1];
        nth_pme_min = -buf_out[2];
        nth_pme_max =  buf_out[3];
    }
    else
#endif
    {
        nth_min     = modth.gnth;
        nth_max     = modth.gnth;
        nth_pme_min = modth.gnth_pme;
        nth_pme_max = modth.gnth_pme;
    }

    /* for group scheme we print PME threads info only */
    if (bFullOmpSupport)
    {
        if (nth_max == nth_min)
        {
            GMX_LOG(mdlog.warning).appendTextFormatted(
                    "Using %d OpenMP thread%s %s",
                    nth_min, nth_min > 1 ? "s" : "",
                    cr->nnodes > 1 ? mpi_str : "");
        }
        else
        {
            GMX_LOG(mdlog.warning).appendTextFormatted(
                    "Using %d - %d OpenMP threads %s",
                    nth_min, nth_max, mpi_str);
        }
    }
    if (bSepPME && (nth_pme_min != nth_min || nth_pme_max != nth_max))
    {
        if (nth_pme_max == nth_pme_min)
        {
            GMX_LOG(mdlog.warning).appendTextFormatted(
                    "Using %d OpenMP thread%s %s for PME",
                    nth_pme_min, nth_pme_min > 1 ? "s" : "",
                    cr->nnodes > 1 ? mpi_str : "");
        }
        else
        {
            GMX_LOG(mdlog.warning).appendTextFormatted(
                    "Using %d - %d OpenMP threads %s for PME",
                    nth_pme_min, nth_pme_max, mpi_str);
        }
    }
    GMX_LOG(mdlog.warning);
}

/*! \brief Detect and warn about oversubscription of cores.
 *
 * \todo This could probably live elsewhere, since it is not specifc
 * to OpenMP, and only needs modth.gnth.
 *
 * \todo Enable this for separate PME nodes as well! */
static void
issueOversubscriptionWarning(const gmx::MDLogger &mdlog,
                             const t_commrec     *cr,
                             int                  nthreads_hw_avail,
                             int                  nppn,
                             gmx_bool             bSepPME)
{
    char sbuf[STRLEN], sbuf1[STRLEN], sbuf2[STRLEN];

    if (bSepPME || 0 != cr->rank_pp_intranode)
    {
        return;
    }

    if (modth.gnth*nppn > nthreads_hw_avail)
    {
        sprintf(sbuf, "threads");
        sbuf1[0] = '\0';
        sprintf(sbuf2, "O");
#if GMX_MPI
        if (modth.gnth == 1)
        {
#if GMX_THREAD_MPI
            sprintf(sbuf, "thread-MPI threads");
#else
            sprintf(sbuf, "MPI processes");
            sprintf(sbuf1, " per rank");
            sprintf(sbuf2, "On rank %d: o", cr->sim_nodeid);
#endif
        }
#endif
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                "WARNING: %sversubscribing the available %d logical CPU cores%s with %d %s.\n"
                "         This will cause considerable performance loss!",
                sbuf2, nthreads_hw_avail, sbuf1, nppn*modth.gnth, sbuf);
    }
}

// Sets the granularity of the nbr. list portions. (Portions are multiples of
// units of size 1/nblDenom.)
const int nblDenom = 1000;

// Mechanism for storing and retrieving computed portions
static int numNbls;
static int numNblsLL;
static std::vector< Range<Ratio> > nblPortionsLocal;
static std::vector< Range<Ratio> > nblPortionsNonlocal;

Range<Ratio> getNblWorkPortion(int nblId, int iloc)
{
    if (iloc == eintLocal) {
        return nblPortionsLocal[nblId];
    }
    else if (iloc == eintNonlocal) {
        return nblPortionsNonlocal[nblId];
    }
    assert(false);
    // Avoid compiler warning
    return {0,1};
}

/*! \brief
 * Return STS overhead for subtask
 *
 * TODO: This is a simple mockup function. Overhead should be timed and
 * recorded for each subtask during execution. For now, we just assume an
 * overhead of 100 microseconds for each nbr. list processed by a subtask,
 * which seems to work rather well.
 */
long getSubTaskOverhead(int st, std::string taskName) {
    const long nblOverhead = 100;
    int firstThreadLL = numNbls - numNblsLL;
    // Currently, the LL threads for local NB have two nbr. lists, while all
    // other threads have one.
    if (st >= firstThreadLL && taskName == "nonbonded_loop_local") {
        return 2*nblOverhead;
    }
    else {
        return nblOverhead;
    }
}

/*! \brief
 * Create a vector of work distributions as units of time (microseconds) for
 * each subtask.
 *
 * The distribution is based on how much time each subtask had available to do
 * work in the previous step.
 *
 * \param[in]  taskName name of task to allocate
 * \param[out] vector to store work distributions
 * \return     total amount of work based on work done in previous step)
 */
long distributeWork(const std::string& taskName, std::vector<long> &stWork) {
    STS* sts = STS::getInstance("force");

    // Get the list of subtasks for the given task
    std::vector<const SubTask*> st = sts->getTask(taskName)->getSubTasks();
    stWork.assign(st.size(),0);

    // Build a sorted vector of subtask start and end times. This divides the
    // run time into intervals, where each interval has a certain subset of
    // threads that are active.
    enum TimeType {start,end};
    struct SubTaskTime {
        TimeType t;
        int i;
        SubTaskTime(TimeType tt, int idx) :t(tt), i(idx) {}
    };
    long remWork = 0;
    std::vector<SubTaskTime> stTimes;
    for (size_t i=0; i<st.size(); i++) {
        stTimes.emplace_back(start,i);
        stTimes.emplace_back(end,i);
        remWork += st[i]->getRunEndTime() - st[i]->getRunStartTime() -
                                            getSubTaskOverhead(i,taskName);
    }
    const long totalWork = remWork;
    std::sort(stTimes.begin(), stTimes.end(), [&st,&taskName](const SubTaskTime& stt1,
                                                              const SubTaskTime& stt2) -> bool
    {
        long t1 = stt1.t == start ?
                              st[stt1.i]->getRunStartTime() : st[stt1.i]->getNextRunAvailTime() -
                                          getSubTaskOverhead(stt1.i,taskName);
        long t2 = stt2.t == start ?
                              st[stt2.i]->getRunStartTime() : st[stt2.i]->getNextRunAvailTime() -
                                          getSubTaskOverhead(stt2.i,taskName);
        return t1 < t2;
    });

    // Iterate through the intervals, "pouring" work into the threads that are
    // active for that interval until all work has been assigned.
    std::set<int> availSTs;

    long istart = 0;
    long iend = st[stTimes[0].i]->getRunStartTime();
    for (size_t i=0; i<2*st.size() - 1; i++) {
        istart = iend;
        iend = stTimes[i+1].t == start ?
               st[stTimes[i+1].i]->getRunStartTime() : st[stTimes[i+1].i]->getNextRunAvailTime() -
                                                       getSubTaskOverhead(stTimes[i+1].i,taskName);
        if (stTimes[i].t == start) {
            availSTs.insert(stTimes[i].i);
        }
        else {
            availSTs.erase(stTimes[i].i);
        }
        long gapWork = std::min((iend - istart)*((long)availSTs.size()), remWork);
        for (int ast : availSTs) {
            stWork[ast] += gapWork / availSTs.size();
        }
        remWork -= gapWork;
        if (remWork <= 0) {
            break;
        }
    }

    return totalWork;
}

//! \brief Helper function to evenly divide a ratio interval among a set of n slots
template<typename Iter>
void dist_work_evenly(Iter start, Iter end, int numerStart, int numerEnd, int denom) {
    int nslots = std::distance(start,end);
    int nunits = numerEnd - numerStart;
    int p = nunits / nslots;
    int r = nunits % nslots;
    int numer = numerStart;
    for (int s=0; s<nslots; s++)
    {
        int inc = s < r ? p+1 : p;
        *(start + s) = { {numer,denom}, {numer+inc,denom} };
        numer += inc;
    }
}

//! \brief Set initial nbr. list portions
void sts_init_nbl_portions() {
    nblPortionsLocal.resize(numNbls);
    nblPortionsNonlocal.resize(numNbls);

    dist_work_evenly(nblPortionsLocal.begin(), nblPortionsLocal.end(), 0, nblDenom, nblDenom);
    dist_work_evenly(nblPortionsNonlocal.begin(), nblPortionsNonlocal.end(), 0, nblDenom, nblDenom);
}

//! \brief Print nbr. list portions
void sts_report_nbl_portions() {
    for (int p = 0; p < numNbls; p++) {
        std::cerr << "Portion " << p << " " << getNblWorkPortion(p, eintLocal).toString() << std::endl;
    }
}

/*! \brief
 * Compute nbr. list portions given a time duration for each subtask
 *
 * \param[in] subTaskTimes vector of time durations in microseconds, one per subtask
 * \param[in] totalWork total amount of work to be done in microseconds
 * \param[in] taskName name of task
 */
void setNblPortions(const std::vector<long> &subTaskTimes, long totalWork, std::string taskName) {
    bool isLocal = taskName == "nonbonded_loop_local" ? true : false;
    std::vector<double> subTaskDist;
    for (long t : subTaskTimes) {
        subTaskDist.push_back(t * 1.0 / totalWork);
    }

    // AU: assigned units. Must assign exactly nblDenom units to nbr. lists
    int totalAUs = 0;
    std::vector<int> au(numNbls,0);

    // First section of subtasks process a single nbr. list
    int numSingles = isLocal ? numNbls-numNblsLL : numNbls;
    for (int i=0; i<numSingles; i++) {
        au[i] = std::max(1.0, round(subTaskDist[i] * nblDenom));
        totalAUs += au[i];
    }

    // Remaining subtasks process two nbr. lists
    for (int i=numSingles,j=numSingles; i<numNbls; i += 2,j++) {
        int units = std::max(1.0, round(subTaskDist[j] * nblDenom / 2));
        au[i] = units;
        au[i+1] = units;
        totalAUs += (au[i] + au[i+1]);
    }

    // Adjust AUs so that exactly nblDenom are assigned
    // Will not terminate if numNbls > nblDenom, due to forced minimum AU of 1
    assert(numNbls <= nblDenom);
    int i=0;
    while (totalAUs != nblDenom) {
        int j = i % numNbls;
        if (totalAUs < nblDenom) {
            au[j]++;
            totalAUs++;
        }
        else if (au[j] > 1) {
            au[j]--;
            totalAUs--;
        }
        i++;
    }

    // Finally, convert AUs to a set of ranges, as expected by the rebalancing code
    std::vector< Range<Ratio> > &v = isLocal ? nblPortionsLocal : nblPortionsNonlocal;
    int numer = 0;
    for (int i=0; i < numNbls; i++) {
        v[i] = { {numer,nblDenom},{numer += au[i],nblDenom} };
    }
}

//! \brief Main entry point to readjust nbr. list sizes (portions)
void sts_adjust_nbl_portions() {
    std::vector<long> subTaskTimes;
    long totalWork = distributeWork("nonbonded_loop_local", subTaskTimes);
    setNblPortions(subTaskTimes, totalWork, "nonbonded_loop_local");
    totalWork = distributeWork("nonbonded_loop_nonlocal", subTaskTimes);
    setNblPortions(subTaskTimes, totalWork, "nonbonded_loop_nonlocal");
}

void sts_report_force_timings() {
    STS* sts = STS::getInstance("force");
    for(int t=0; t<STS::getNumThreads(); t++) {
        dprintf(2, "Thread %d\n", t);
        for (int s = 0; s<sts->getNumSubTasks(t); s++) {
            const SubTask* st = sts->getSubTask(t,s);
            std::cerr << st->getTask()->getLabel() << " " << st->getWaitStartTime() << " " << st->getRunStartTime() << " " << st->getRunEndTime() << std::endl;
            if (t==0 && st->getTask()->getLabel() == "PME1") {
                int comm_id = 0;
                for (long t : st->getAuxTimes("fft_comm")) {
                    std::cerr << "fft_comm" << comm_id << " " << t << std::endl;
                    comm_id++;
                }
            }
        }
    }
}

//! \brief Creates the STS schedule for computing forces.
void sts_create_schedule_for_force_compute(const t_commrec *cr, bool doInit, SchedType stype, double ratioThreadsForNB)
{
    // Initialize if and only if this is the first call. Note that doInit is
    // not really needed but makes the difference between the first call and
    // subsequent calls explicit.
    static bool firstCall = true;
    assert(firstCall == doInit);
    firstCall = false;

    // Ratio and schedule type are set permanently by the first call.
    // Can only change thread counts before thread initialization (and before MD starts)
    // Thread counts for nonbonded, bonded, and PME cannot change during run.
    static SchedType schedType = SchedType::DEFAULT;
    if (doInit) {
        schedType = stype;
    }

    // Schedule for force computations. Note that STS caches the instance for us.
    STS* sts;
    if (doInit) {
        sts = new STS("force");
    }
    else {
        sts = STS::getInstance("force");
    }
    sts->clearAssignments();

    const int nthreads = STS::getNumThreads();
    if (schedType == SchedType::DEFAULT || nthreads == 1) {
        sts->setDefaultSchedule();
        if (doInit) {
            new MMBarrier(nthreads, "fft");
        }
        // TODO: We should also set the env. variables in this case, even
        // though the current default is to use all threads.
        return;
    }

    // Various constants

    // Threads that do non-NB force computations (low-level (LL) forces)
    const int  mainThreadLL      = 0;
    const int  nthreadsLL        = std::max(2.0,(1.0-ratioThreadsForNB) * nthreads);
    // Note that main NB thread is not thread 1 but first non-LL thread
    const int  mainThreadNB     = nthreadsLL;
    // All threads except 0 do some NB work
    const int  nthreadsNB = nthreads-1;
    const int nthreadsOnlyNB = nthreads-nthreadsLL;
    // Threads that do both LL and NB (threadsLL without thread 0)
    const int  nthreadsBoth = nthreadsLL-1;
    const bool doPP  = (cr->duty & DUTY_PP);
    const bool doPME = (cr->duty & DUTY_PME);
    // PME is done by the lowlevel threads unless this is a PME-only process.
    // Then all threads are used for PME.
    const int firstThreadPME = doPP ? mainThreadLL : 0;
    const int nthreadsPME = doPME ? (doPP ? nthreadsLL : nthreads) : 0;
    // One pair list per NB thread, except for "both" threads to support
    // alternating between NB and LL work (could be avoided if coroutines
    // were supported). These threads receive three additional pair lists.
    const int nPairLists = nthreadsNB + 3*nthreadsBoth;

    numNbls = nPairLists;
    numNblsLL = 4*nthreadsBoth;
    sts_init_nbl_portions();

    // Various thread sets

    std::vector<int> threadsLL(nthreadsLL);
    std::iota(threadsLL.begin(),threadsLL.end(), mainThreadLL);
    std::vector<int> threadsNB(nthreadsNB);
    std::iota(threadsNB.begin(), threadsNB.end(), 1);
    std::vector<int> threadsOnlyNB(nthreadsOnlyNB);
    std::iota(threadsOnlyNB.begin(), threadsOnlyNB.end(), mainThreadNB);
    std::vector<int> threadsBoth(nthreadsBoth);
    std::iota(threadsBoth.begin(), threadsBoth.end(), mainThreadLL+1);
    std::vector<int> threadsPME(nthreadsPME);
    std::iota(threadsPME.begin(), threadsPME.end(), firstThreadPME);

    if (doPP) {
        if (DOMAINDECOMP(cr)) {
            sts->assign_run("xcomm", mainThreadLL);
            sts->assign_loop("nbat_x_to_x", mainThreadLL, {0,1});
        }

        // Assign NB-only threads. Also, "both" threads start on local NB and
        // so are assigned here.
        sts->assign_run("nonbonded", mainThreadNB);
        sts->assign_loop("nonbonded_loop_local", threadsOnlyNB,
                                                 { 0, {nthreadsOnlyNB,nPairLists} });
        sts->assign_loop("nonbonded_loop_local", threadsBoth,
                                                 { {nthreadsOnlyNB,nPairLists},
                                                   {nthreadsOnlyNB+2*nthreadsBoth,nPairLists} });
        if (DOMAINDECOMP(cr)) {
            sts->assign_loop("nonbonded_loop_nonlocal", threadsOnlyNB,
                                                        { 0, {nthreadsOnlyNB,nPairLists} });
        }
        if (nPairLists > 1) {
            sts->assign_loop("stdreduce", threadsOnlyNB);
            // sts->assign_loop("nbat_f_to_f", threadsNB);
        }
        
        sts->assign_run("lowlevel", 0);
        sts->assign_loop("listed_forces2", threadsLL);
        sts->assign_loop("listed_forces1", threadsLL);
        // Local NB during redistribution prior to PME
        Range<Ratio> nbRange = { {nthreadsOnlyNB+2*nthreadsBoth,nPairLists},1 };
        sts->assign_loop("nonbonded_loop_local", threadsBoth, nbRange);
    }

    if (doPME) {
        // If PP done, PME will be run as part of the lowlevel forces.
        if (!doPP) {
            sts->assign_run("PME", 0);
        }
        sts->assign_loop("calc_splines", threadsPME);
        sts->assign_loop("pme_spread", threadsPME);
        if (nthreadsPME > 1) {
            sts->assign_loop("pme_spread2", threadsPME);
        }

        sts->assign_loop("PME1", threadsPME);
        if (doPP) {
            // Assign nonlocal NB iterations to "both" threads. This work will be done
            // during PME communications.
            for (int numer = nthreadsOnlyNB; numer < nPairLists; numer += nthreadsBoth) {
                Range<Ratio> nbRange = { {numer, nPairLists},
                                         {numer + nthreadsBoth, nPairLists} };
                sts->assign_loop("nonbonded_loop_nonlocal", threadsBoth, nbRange);
            }
            sts->setHighPriority("nonbonded_loop_nonlocal");
        }
        sts->assign_loop("gather_f_bsplines", threadsPME);
        // sts->assign_loop("PME2", threadsPME);
        // sts->assign_loop("PME3", threadsPME);
        // sts->assign_loop("PME4", threadsPME);
        // sts->assign_loop("corr", threadsPME);

        // Create STS barrier used by FFT library.
        // We create it here because we know the number of threads.
        if (doInit) {
            new MMBarrier(nthreadsPME, "fft");
        }
    }

    // Tasks not yet assigned (done either serially or with all threads if the
    // default schedule is used).
    // listed_manage1, listed_manage2, nbnxn_grid
    // nbnxn_search1, nbnxn_search2, nbnxn_search3, nbnxn_search4, nbnxn_search5
    // calc_splines, pme_spread, pme_spread2, gather_f_bsplines

    if (!doInit) {
        return;
    }

    // Adjust thread counts to match STS schedule
    // For portability, use putenv and use static buffers that are never deallocated.
    // (See putenv man page for details)
    static char nb_env_string[40];
    strcpy(nb_env_string, "GMX_NONBONDED_NUM_THREADS=");
    strcat(nb_env_string, std::to_string(nPairLists).c_str());
    putenv(nb_env_string);

    static char listed_env_string[40];
    strcpy(listed_env_string, "GMX_LISTED_FORCES_NUM_THREADS=");
    strcat(listed_env_string, std::to_string(nthreadsLL).c_str());
    putenv(listed_env_string);

    static char pme_env_string[40];
    strcpy(pme_env_string, "GMX_PME_NUM_THREADS=");
    strcat(pme_env_string, std::to_string(nthreadsPME).c_str());
    putenv(pme_env_string);
}

void gmx_omp_nthreads_init(const gmx::MDLogger &mdlog, t_commrec *cr,
                           int nthreads_hw_avail,
                           int omp_nthreads_req,
                           int omp_nthreads_pme_req,
                           gmx_bool bThisNodePMEOnly,
                           gmx_bool bFullOmpSupport)
{
    int        nppn;
    gmx_bool   bSepPME;

    const bool bOMP = GMX_STS;

    /* number of MPI processes/threads per physical node */
    nppn = cr->nrank_intranode;

    bSepPME = ( (cr->duty & DUTY_PP) && !(cr->duty & DUTY_PME)) ||
        (!(cr->duty & DUTY_PP) &&  (cr->duty & DUTY_PME));

    manage_number_of_openmp_threads(mdlog, cr, bOMP,
                                    nthreads_hw_avail,
                                    omp_nthreads_req, omp_nthreads_pme_req,
                                    bThisNodePMEOnly, bFullOmpSupport,
                                    nppn, bSepPME);
#if GMX_THREAD_MPI
    /* Non-master threads have to wait for the OpenMP management to be
     * done, so that code elsewhere that uses OpenMP can be certain
     * the setup is complete. */
    if (PAR(cr))
    {
        MPI_Barrier(cr->mpi_comm_mysim);
    }
#endif

    reportOpenmpSettings(mdlog, cr, bOMP, bFullOmpSupport, bSepPME);
    issueOversubscriptionWarning(mdlog, cr, nthreads_hw_avail, nppn, bSepPME);
}

int gmx_omp_nthreads_get(int mod)
{
    if (mod < 0 || mod >= emntNR)
    {
        /* invalid module queried */
        return -1;
    }
    else
    {
        return modth.nth[mod];
    }
}

void
gmx_omp_nthreads_set(int mod, int nthreads)
{
    /* Catch an attempt to set the number of threads on an invalid
     * OpenMP module. */
    GMX_RELEASE_ASSERT(mod >= 0 && mod < emntNR, "Trying to set nthreads on invalid OpenMP module");

    modth.nth[mod] = nthreads;
}
