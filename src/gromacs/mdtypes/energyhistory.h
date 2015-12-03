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
#ifndef GMX_MDLIB_ENERGYHISTORY_H
#define GMX_MDLIB_ENERGYHISTORY_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

/* energy history for delta_h histograms */
typedef struct delta_h_history_t
{
    int      nndh;             /* the number of energy difference lists */
    int     *ndh;              /* the number in each energy difference list */
    real   **dh;               /* the energy difference lists */

    double   start_time;       /* the start time of these energy diff blocks */
    double   start_lambda;     /* lambda at start time */

    gmx_bool start_lambda_set; /* whether the lambda value is set. Here
                                  For backward-compatibility. */
} delta_h_history_t;


typedef struct energyhistory_t
{
    gmx_int64_t        nsteps;       /* The number of steps in the history            */
    gmx_int64_t        nsum;         /* The nr. of steps in the ener_ave and ener_sum */
    double         *   ener_ave;     /* Energy term history sum to get fluctuations   */
    double         *   ener_sum;     /* Energy term history sum to get fluctuations   */
    int                nener;        /* Number of energy terms in two previous arrays */
    gmx_int64_t        nsteps_sim;   /* The number of steps in ener_sum_sim      */
    gmx_int64_t        nsum_sim;     /* The number of frames in ener_sum_sim     */
    double         *   ener_sum_sim; /* Energy term history sum of the whole sim      */

    delta_h_history_t *dht;          /* The BAR energy differences */
}
energyhistory_t;

/* \brief Initialize an energy history structure
 */
void init_energyhistory(energyhistory_t * enerhist);

/* \brief Destroy an energy history structure
 */
void done_energyhistory(energyhistory_t * enerhist);

#endif
