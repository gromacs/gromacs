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

#ifndef GMX_OMP_NTHREADS_H
#define GMX_OMP_NTHREADS_H

#include <cstdio>

#include "gromacs/utility/basedefinitions.h"

struct t_commrec;

namespace gmx
{
class MDLogger;
}

/** Enum values corresponding to multithreaded algorithmic modules. */
enum class ModuleMultiThread : int
{
    /* Default is meant to be used in OMP regions outside the named
     * algorithmic modules listed below. */
    Default,
    Domdec,
    Pairsearch,
    Nonbonded,
    Bonded,
    Pme,
    Update,
    VirtualSite,
    Lincs,
    Settle,
    Count
};

/*! \brief
 * Initializes the per-module thread count.
 *
 * It is compatible with tMPI, thread-safety is ensured (for the features
 * available with tMPI).
 * This function should be called only once during the initialization of mdrun. */
void gmx_omp_nthreads_init(const gmx::MDLogger& fplog,
                           t_commrec*           cr,
                           int                  maxThreads,
                           int                  numRanksOnThisNode,
                           int                  omp_nthreads_req,
                           int                  omp_nthreads_pme_req,
                           gmx_bool             bCurrNodePMEOnly);

/*! \brief
 * Returns the number of threads to be used in the given module \p mod. */
int gmx_omp_nthreads_get(ModuleMultiThread mod);

/*! \brief
 * Returns the number of threads to be used in the given module \p mod for simple rvec operations.
 *
 * When the, potentially, parallel task only consists of a loop of clear_rvec
 * or rvec_inc for nrvec elements, the OpenMP overhead might be higher than
 * the reduction in computational cost due to parallelization. This routine
 * returns 1 when the overhead is expected to be higher than the gain.
 */
static inline int gmx_omp_nthreads_get_simple_rvec_task(ModuleMultiThread mod, int nrvec)
{
    /* There can be a relatively large overhead to an OpenMP parallel for loop.
     * This overhead increases, slowly, with the numbe of threads used.
     * The computational gain goes as 1/#threads. The two effects combined
     * lead to a cross-over point for a (non-)parallel loop at loop count
     * that is not strongly dependent on the thread count.
     * Note that a (non-)parallel loop can have benefit later in the code
     * due to generating more cache hits, depending on how the next lask
     * that accesses the same data is (not) parallelized over threads.
     *
     * A value of 2000 is the switch-over point for Haswell without
     * hyper-threading. With hyper-threading it is about a factor 1.5 higher.
     */
    const int nrvec_omp = 2000;

    if (nrvec < nrvec_omp)
    {
        return 1;
    }
    else
    {
        return gmx_omp_nthreads_get(mod);
    }
}

/*! \brief Sets the number of threads to be used in module.
 *
 * Intended for use in testing. */
void gmx_omp_nthreads_set(ModuleMultiThread mod, int nthreads);

/*! \brief
 * Read the OMP_NUM_THREADS env. var. and check against the value set on the
 * command line. */
void gmx_omp_nthreads_read_env(const gmx::MDLogger& mdlog, int* nthreads_omp);

#endif
