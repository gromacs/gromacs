/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef GMX_OMP_NTHREADS
#define GMX_OMP_NTHREADS
#include "visibility.h"
#include "types/commrec.h"


/*! Enum values corresponding to multithreaded algorithmic modules. */
typedef enum module_nth
{
    /* Default is meant to be used in OMP regions outside the named
     * algorithmic modules listed below. */
    emntDefault, emntDomdec, emntPairsearch, emntNonbonded,
    emntBonded, emntPME,  emntUpdate, emntVSITE, emntLINCS, emntSETTLE,
    emntNR
} module_nth_t;

/*! Initializes the per-module thread count. It is compatible with tMPI,
 *  thread-safety is ensured (for the features available with tMPI).
 *  This function should caled only once during the initialization of mdrun. */
GMX_LIBGMX_EXPORT
void gmx_omp_nthreads_init(FILE *fplog, t_commrec *cr,
                           int nthreads_hw_avail,
                           int omp_nthreads_req,
                           int omp_nthreads_pme_req,
                           gmx_bool bCurrNodePMEOnly,
                           gmx_bool bFullOmpSupport);

/*! Returns the number of threads to be used in the given module m. */
GMX_LIBGMX_EXPORT
int gmx_omp_nthreads_get(int mod);

/*! Read the OMP_NUM_THREADS env. var. and check against the value set on the command line. */
GMX_LIBGMX_EXPORT
void gmx_omp_nthreads_read_env(int     *nthreads_omp,
                               gmx_bool bIsSimMaster);

#endif /* GMX_OMP_NTHREADS */
