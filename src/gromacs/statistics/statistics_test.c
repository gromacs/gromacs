/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2011,2014, by the GROMACS development team, led by
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

#include <stdio.h>

#include "gromacs/math/vec.h"
#include "gromacs/random/random.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

static void horizontal()
{
    gmx_rng_t   rng;
    gmx_stats_t straight;
    int         i, ok, n = 1000;
    real        y, a, b, da, db, aver, sigma, error, chi2, R, *xh, *yh;
    FILE       *fp;

    rng      = gmx_rng_init(13);
    straight = gmx_stats_init();
    for (i = 0; (i < n); i++)
    {
        y = gmx_rng_uniform_real(rng);
        if ((ok = gmx_stats_add_point(straight, i, y, 0, 0)) != estatsOK)
        {
            fprintf(stderr, "%s\n", gmx_stats_message(ok));
        }
    }
    /* Horizontal test */
    if ((ok = gmx_stats_get_ase(straight, &aver, &sigma, &error)) != estatsOK)
    {
        fprintf(stderr, "%s\n", gmx_stats_message(ok));
    }
    fp = fopen("straight.xvg", "w");
    if ((ok = gmx_stats_dump_xy(straight, fp)) != estatsOK)
    {
        fprintf(stderr, "%s\n", gmx_stats_message(ok));
    }
    fclose(fp);
    printf("Horizontal line: average %g, sigma %g, error %g\n", aver, sigma, error);
    if ((ok = gmx_stats_done(straight)) != estatsOK)
    {
        fprintf(stderr, "%s\n", gmx_stats_message(ok));
    }
}

static void line()
{
    gmx_rng_t   rng;
    gmx_stats_t line;
    int         i, dy, ok, n = 1000;
    real        y, a, b, da, db, aver, sigma, error, chi2, R, rfit;
    const real  a0 = 0.23, b0 = 2.7;
    FILE       *fp;

    for (dy = 0; (dy < 2); dy++)
    {
        rng      = gmx_rng_init(13);
        line     = gmx_stats_init();
        for (i = 0; (i < n); i++)
        {
            y = a0*i+b0+50*(gmx_rng_uniform_real(rng)-0.5);
            if ((ok = gmx_stats_add_point(line, i, y, 0, dy*0.1)) != estatsOK)
            {
                fprintf(stderr, "%s\n", gmx_stats_message(ok));
            }
        }
        /* Line with slope test */
        if ((ok = gmx_stats_get_ab(line, elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &rfit)) != estatsOK)
        {
            fprintf(stderr, "%s\n", gmx_stats_message(ok));
        }
        if ((ok = gmx_stats_get_corr_coeff(line, &R)) != estatsOK)
        {
            fprintf(stderr, "%s\n", gmx_stats_message(ok));
        }
        if (dy == 0)
        {
            fp = fopen("line0.xvg", "w");
        }
        else
        {
            fp = fopen("line1.xvg", "w");
        }
        if ((ok = gmx_stats_dump_xy(line, fp)) != estatsOK)
        {
            fprintf(stderr, "%s\n", gmx_stats_message(ok));
        }
        fclose(fp);
        printf("Line with eqn. y = %gx + %g with noise%s\n", a0, b0,
               (dy == 0) ? "" : " and uncertainties");
        printf("Found: a = %g +/- %g, b = %g +/- %g\n", a, da, b, db);
        if ((ok = gmx_stats_done(line)) != estatsOK)
        {
            fprintf(stderr, "%s\n", gmx_stats_message(ok));
        }
        gmx_rng_destroy(rng);
    }
}

static void histogram()
{
    gmx_rng_t   rng;
    gmx_stats_t camel;
    int         i, ok, n = 1000, norm;
    real        y, a, b, da, db, aver, sigma, error, chi2, R, *xh, *yh;
    const real  a0 = 0.23, b0 = 2.7;
    FILE       *fp;
    char        fn[256];

    for (norm = 0; (norm < 2); norm++)
    {
        rng      = gmx_rng_init(13);
        camel    = gmx_stats_init();
        for (i = 0; (i < n); i++)
        {
            y = sqr(gmx_rng_uniform_real(rng));
            if ((ok = gmx_stats_add_point(camel, i, y+1, 0, 0)) != estatsOK)
            {
                fprintf(stderr, "%s\n", gmx_stats_message(ok));
            }
            y = sqr(gmx_rng_uniform_real(rng));
            if ((ok = gmx_stats_add_point(camel, i+0.5, y+2, 0, 0)) != estatsOK)
            {
                fprintf(stderr, "%s\n", gmx_stats_message(ok));
            }
        }
        /* Histogram test */
        if ((ok = gmx_stats_make_histogram(camel, 0, 101, norm, &xh, &yh)) != estatsOK)
        {
            fprintf(stderr, "%s\n", gmx_stats_message(ok));
        }
        sprintf(fn, "histo%d-data.xvg", norm);
        fp = fopen(fn, "w");
        gmx_stats_dump_xy(camel, fp);
        fclose(fp);
        sprintf(fn, "histo%d.xvg", norm);
        fp = fopen(fn, "w");
        for (i = 0; (i < 101); i++)
        {
            fprintf(fp, "%12g  %12g\n", xh[i], yh[i]);
        }
        fclose(fp);
        sfree(xh);
        sfree(yh);
    }
}

int main(int argc, char *argv[])
{
    line();
    horizontal();
    histogram();

    return 0;
}
