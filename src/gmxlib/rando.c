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

#include <time.h>
#ifdef GMX_NATIVE_WINDOWS
#include <process.h>
#endif
#include "sysstuff.h"
#include "typedefs.h"
#include "random.h"

int make_seed(void)
{
#ifdef GMX_NATIVE_WINDOWS
    return (int)_getpid();
#else
    return (int)getpid();
#endif
}

real rando(int *ig)
/* generate a random number. */
{
    int  irand;

    int  m    = 100000000;
    real rm   = 100000000.0; /* same number as m, but real format */
    int  m1   = 10000;
    int  mult = 31415821;

    real r;
    int  irandh, irandl, multh, multl;

    irand = abs(*ig) % m;

    /* multiply irand by mult, but take into account that overflow
     * must be discarded, and do not generate an error.
     */
    irandh = irand / m1;
    irandl = irand % m1;
    multh  = mult / m1;
    multl  = mult % m1;
    irand  = ((irandh*multl+irandl*multh) % m1) * m1 + irandl*multl;
    irand  = (irand + 1) % m;

    /* convert irand to a real random number between 0 and 1. */
    r = (irand / 10);
    r = r * 10 / rm;
    if ((r <= 0) || (r > 1))
    {
        r = 0.0;
    }
    *ig = irand;

    return r;
}
