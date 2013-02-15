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

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "vec.h"
#include "coulomb.h"
#include "smalloc.h"
#include "physics.h"
#include "txtdump.h"
#include "futil.h"
#include "names.h"
#include "writeps.h"
#include "macros.h"
#include "xvgr.h"
#include "gmxfio.h"

#ifdef GMX_THREAD_MPI
#include "thread_mpi/threads.h"
#endif

#define p2(x) ((x)*(x))
#define p3(x) ((x)*(x)*(x))
#define p4(x) ((x)*(x)*(x)*(x))

static real                A, A_3, B, B_4, C, c1, c2, c3, c4, c5, c6, One_4pi, FourPi_V, Vol, N0;
#ifdef GMX_THREAD_MPI
static tMPI_Thread_mutex_t shift_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
#endif


void set_shift_consts(FILE *log, real r1, real rc, rvec box, t_forcerec *fr)
{
#ifdef GMX_THREAD_MPI
    /* at the very least we shouldn't allow multiple threads to set these
       simulataneously */
    tMPI_Thread_mutex_lock(&shift_mutex);
#endif
    /* A, B and C are recalculated in tables.c */
    if (r1 < rc)
    {
        A   = (2*r1-5*rc)/(p3(rc)*p2(rc-r1));
        B   = (4*rc-2*r1)/(p3(rc)*p3(rc-r1));
        /*C   = (10*rc*rc-5*rc*r1+r1*r1)/(6*rc*rc); Hermans Eq. not correct */
    }
    else
    {
        gmx_fatal(FARGS, "r1 (%f) >= rc (%f) in %s, line %d",
                  r1, rc, __FILE__, __LINE__);
    }

    A_3 = A/3.0;
    B_4 = B/4.0;
    C   = 1/rc-A_3*p3(rc-r1)-B_4*p4(rc-r1);
    N0  = 2.0*M_PI*p3(rc)*p3(rc-r1);

    Vol      = (box[XX]*box[YY]*box[ZZ]);
    FourPi_V = 4.0*M_PI/Vol;

    if (debug)
    {
        fprintf(debug, "Constants for short-range and fourier stuff:\n"
                "r1 = %10.3f,  rc = %10.3f\n"
                "A  = %10.3e,  B  = %10.3e,  C  = %10.3e, FourPi_V = %10.3e\n",
                r1, rc, A, B, C, FourPi_V);
    }

    /* Constants derived by Mathematica */
    c1 = -40*rc*rc    + 50*rc*r1    - 16*r1*r1;
    c2 =  60*rc       - 30*r1;
    c3 = -10*rc*rc*rc + 20*rc*rc*r1 - 13*rc*r1*r1 + 3*r1*r1*r1;
    c4 = -20*rc*rc    + 40*rc*r1    - 14*r1*r1;
    c5 = -c2;
    c6 = -5*rc*rc*r1  +  7*rc*r1*r1 - 2*r1*r1*r1;

    if (debug)
    {
        fprintf(debug, "c1 = %10.3e,  c2 = %10.3e,  c3 = %10.3e\n"
                "c4 = %10.3e,  c5 = %10.3e,  c6 = %10.3e,  N0 = %10.3e\n",
                c1, c2, c3, c4, c5, c6, N0);
    }

    One_4pi = 1.0/(4.0*M_PI);
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_unlock(&shift_mutex);
#endif
}
