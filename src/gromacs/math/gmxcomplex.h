/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
#ifndef GMX_MATH_GMXCOMPLEX_H
#define GMX_MATH_GMXCOMPLEX_H

#include <math.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

typedef struct {
    real re, im;
} t_complex;

typedef t_complex cvec[DIM];

static t_complex rcmul(real r, t_complex c)
{
    t_complex d;

    d.re = r*c.re;
    d.im = r*c.im;

    return d;
}

static t_complex rcexp(real r)
{
    t_complex c;

    c.re = (real)cos(r);
    c.im = (real)sin(r);

    return c;
}


static t_complex cadd(t_complex a, t_complex b)
{
    t_complex c;

    c.re = a.re+b.re;
    c.im = a.im+b.im;

    return c;
}

static t_complex csub(t_complex a, t_complex b)
{
    t_complex c;

    c.re = a.re-b.re;
    c.im = a.im-b.im;

    return c;
}

static t_complex cmul(t_complex a, t_complex b)
{
    t_complex c;

    c.re = a.re*b.re - a.im*b.im;
    c.im = a.re*b.im + a.im*b.re;

    return c;
}

static t_complex conjugate(t_complex c)
{
    t_complex d;

    d.re =  c.re;
    d.im = -c.im;

    return d;
}

static real cabs2(t_complex c)
{
    real abs2;
    abs2 = (c.re*c.re)+(c.im*c.im);

    return abs2;
}



static t_complex cdiv(t_complex teller, t_complex noemer)
{
    t_complex res, anoemer;

    anoemer = cmul(conjugate(noemer), noemer);
    res     = cmul(teller, conjugate(noemer));

    return rcmul(1.0/anoemer.re, res);
}
#endif
