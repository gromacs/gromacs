/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "statutil.h"
#include "gmx_fatal.h"

#ifdef GMX_THREAD_MPI
#include "thread_mpi.h"
#endif

/* The source code in this file should be thread-safe.
         Please keep it that way. */

/* Globals for trajectory input */
typedef struct {
    real     t;
    gmx_bool bSet;
} t_timecontrol;

static t_timecontrol timecontrol[TNR] = {
    { 0, FALSE },
    { 0, FALSE },
    { 0, FALSE }
};

#ifdef GMX_THREAD_MPI
static tMPI_Thread_mutex_t tc_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
#endif

gmx_bool bTimeSet(int tcontrol)
{
    gmx_bool ret;

#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&tc_mutex);
#endif
    range_check(tcontrol, 0, TNR);
    ret = timecontrol[tcontrol].bSet;
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_unlock(&tc_mutex);
#endif

    return ret;
}

real rTimeValue(int tcontrol)
{
    real ret;

#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&tc_mutex);
#endif
    range_check(tcontrol, 0, TNR);
    ret = timecontrol[tcontrol].t;
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_unlock(&tc_mutex);
#endif
    return ret;
}

void setTimeValue(int tcontrol, real value)
{
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&tc_mutex);
#endif
    range_check(tcontrol, 0, TNR);
    timecontrol[tcontrol].t    = value;
    timecontrol[tcontrol].bSet = TRUE;
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_unlock(&tc_mutex);
#endif
}
