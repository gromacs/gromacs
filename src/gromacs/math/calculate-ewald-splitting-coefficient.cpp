/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "calculate-ewald-splitting-coefficient.h"

#include <cmath>

#include "gromacs/math/utilities.h"
#include "gromacs/utility/real.h"

real calc_ewaldcoeff_q(real rc, real rtol)
{
    real beta = 5, low, high;
    int  n, i = 0;

    do
    {
        i++;
        beta *= 2;
    }
    while (gmx_erfc(beta*rc) > rtol);

    /* Do a binary search with tolerance 2^-60 */
    n    = i+60;
    low  = 0;
    high = beta;
    for (i = 0; i < n; i++)
    {
        beta = (low+high)/2;
        if (gmx_erfc(beta*rc) > rtol)
        {
            low = beta;
        }
        else
        {
            high = beta;
        }
    }
    return beta;
}

static real compute_lj_function(real beta, real rc)
{
    real xrc, xrc2, xrc4, result;
    xrc    = beta*rc;
    xrc2   = xrc*xrc;
    xrc4   = xrc2*xrc2;
    result = std::exp(-xrc2)*(1 + xrc2 + xrc4/2.0);

    return result;
}

real calc_ewaldcoeff_lj(real rc, real rtol)
{
    real beta = 5, low, high;
    int  n, i = 0;

    do
    {
        i++;
        beta *= 2.0;
    }
    while (compute_lj_function(beta, rc) > rtol);

    /* Do a binary search with tolerance 2^-60 */
    n    = i + 60;
    low  = 0;
    high = beta;
    for (i = 0; i < n; ++i)
    {
        beta = (low + high) / 2.0;
        if (compute_lj_function(beta, rc) > rtol)
        {
            low = beta;
        }
        else
        {
            high = beta;
        }
    }
    return beta;
}
