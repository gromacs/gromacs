/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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
/*! \brief
 * Implement routines for integrating a data set
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "integrate.h"

#include <stdio.h>

#include <cmath>

#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"

/*! \brief Integrate a function and printe the integral value. */
real print_and_integrate(FILE *fp, int n, real dt, const real c[],
                         const real *fit, int nskip)
{
    real c0, sum;
    int  j;

    /* Use trapezoidal rule for calculating integral */
    sum = 0.0;
    for (j = 0; (j < n); j++)
    {
        c0 = c[j];
        if (fp && (nskip == 0 || j % nskip == 0))
        {
            fprintf(fp, "%10.3f  %10.5f\n", j*dt, c0);
        }
        if (j > 0)
        {
            sum += dt*(c0+c[j-1]);
        }
    }
    if (fp)
    {
        fprintf(fp, "&\n");
        if (fit)
        {
            for (j = 0; (j < n); j++)
            {
                if (nskip == 0 || j % nskip == 0)
                {
                    fprintf(fp, "%10.3f  %10.5f\n", j*dt, fit[j]);
                }
            }
            fprintf(fp, "&\n");
        }
    }
    return sum*0.5;
}

/*! \brief Compute and return the integral of a function. */
real evaluate_integral(int n, const real x[], const real y[],
                       const real dy[], real aver_start,
                       real *stddev)
{
    double sum, sum_var, w;
    double sum_tail = 0, sum2_tail = 0;
    int    j, nsum_tail = 0;

    /* Use trapezoidal rule for calculating integral */
    if (n <= 0)
    {
        gmx_fatal(FARGS, "Evaluating integral: n = %d (file %s, line %d)",
                  n, __FILE__, __LINE__);
    }

    sum     = 0;
    sum_var = 0;
    for (j = 0; (j < n); j++)
    {
        w = 0;
        if (j > 0)
        {
            w += 0.5*(x[j] - x[j-1]);
        }
        if (j < n-1)
        {
            w += 0.5*(x[j+1] - x[j]);
        }
        sum += w*y[j];
        if (dy)
        {
            /* Assume all errors are uncorrelated */
            sum_var += gmx::square(w*dy[j]);
        }

        if ((aver_start > 0) && (x[j] >= aver_start))
        {
            sum_tail  += sum;
            sum2_tail += std::sqrt(sum_var);
            nsum_tail += 1;
        }
    }

    if (nsum_tail > 0)
    {
        sum = sum_tail/nsum_tail;
        /* This is a worst case estimate, assuming all stddev's are correlated. */
        *stddev = sum2_tail/nsum_tail;
    }
    else
    {
        *stddev = std::sqrt(sum_var);
    }

    return sum;
}
