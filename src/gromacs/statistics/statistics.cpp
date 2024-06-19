/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "statistics.h"

#include <cmath>
#include <cstdlib>

#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"


typedef struct gmx_stats
{
    double  aa, a, b, sigma_a, sigma_b, aver, sigma_aver, error;
    double  rmsd, Rfit, chi2;
    double *x, *y, *dx, *dy;
    int     computed;
    int     np, nalloc;
} gmx_stats;

gmx_stats_t gmx_stats_init()
{
    gmx_stats* stats;

    snew(stats, 1);

    return static_cast<gmx_stats_t>(stats);
}

void gmx_stats_free(gmx_stats_t gstats)
{
    gmx_stats* stats = static_cast<gmx_stats*>(gstats);

    sfree(stats->x);
    sfree(stats->y);
    sfree(stats->dx);
    sfree(stats->dy);
    sfree(stats);
}

void gmx_stats_add_point(gmx_stats_t gstats, double x, double y, double dx, double dy)
{
    gmx_stats* stats = gstats;

    if (stats->np + 1 >= stats->nalloc)
    {
        if (stats->nalloc == 0)
        {
            stats->nalloc = 1024;
        }
        else
        {
            stats->nalloc *= 2;
        }
        srenew(stats->x, stats->nalloc);
        srenew(stats->y, stats->nalloc);
        srenew(stats->dx, stats->nalloc);
        srenew(stats->dy, stats->nalloc);
        for (int i = stats->np; (i < stats->nalloc); i++)
        {
            stats->x[i]  = 0;
            stats->y[i]  = 0;
            stats->dx[i] = 0;
            stats->dy[i] = 0;
        }
    }
    stats->x[stats->np]  = x;
    stats->y[stats->np]  = y;
    stats->dx[stats->np] = dx;
    stats->dy[stats->np] = dy;
    stats->np++;
    stats->computed = 0;
}

static void gmx_stats_compute(gmx_stats* stats, int weight)
{
    double yx, xx, sx, sy, dy, chi2, d2;
    double ssxx, ssyy, ssxy;
    double w, wtot, yx_nw, sy_nw, sx_nw, yy_nw, xx_nw, dx2;

    int N = stats->np;

    if (stats->computed == 0)
    {
        GMX_RELEASE_ASSERT(N >= 1, "Must have points to work on");

        xx = xx_nw = 0;
        yy_nw      = 0;
        yx = yx_nw = 0;
        sx = sx_nw = 0;
        sy = sy_nw = 0;
        wtot       = 0;
        d2         = 0;
        for (int i = 0; (i < N); i++)
        {
            d2 += gmx::square(stats->x[i] - stats->y[i]);
            if (((stats->dy[i]) != 0.0) && (weight == elsqWEIGHT_Y))
            {
                w = 1 / gmx::square(stats->dy[i]);
            }
            else
            {
                w = 1;
            }

            wtot += w;

            xx += w * gmx::square(stats->x[i]);
            xx_nw += gmx::square(stats->x[i]);

            yy_nw += gmx::square(stats->y[i]);

            yx += w * stats->y[i] * stats->x[i];
            yx_nw += stats->y[i] * stats->x[i];

            sx += w * stats->x[i];
            sx_nw += stats->x[i];

            sy += w * stats->y[i];
            sy_nw += stats->y[i];
        }

        /* Compute average, sigma and error */
        stats->aver       = sy_nw / N;
        stats->sigma_aver = std::sqrt(yy_nw / N - gmx::square(sy_nw / N));
        stats->error      = stats->sigma_aver / std::sqrt(static_cast<double>(N));

        /* Compute RMSD between x and y */
        stats->rmsd = std::sqrt(d2 / N);

        /* Correlation coefficient for data */
        yx_nw /= N;
        xx_nw /= N;
        yy_nw /= N;
        sx_nw /= N;
        sy_nw /= N;
        ssxx = N * (xx_nw - gmx::square(sx_nw));
        ssyy = N * (yy_nw - gmx::square(sy_nw));
        ssxy = N * (yx_nw - (sx_nw * sy_nw));

        /* Compute straight line through datapoints, either with intercept
           zero (result in aa) or with intercept variable (results in a
           and b) */
        yx = yx / wtot;
        xx = xx / wtot;
        sx = sx / wtot;
        sy = sy / wtot;

        stats->aa = (yx / xx);
        stats->a  = (yx - sx * sy) / (xx - sx * sx);
        stats->b  = (sy) - (stats->a) * (sx);

        /* Compute chi2, deviation from a line y = ax+b. */
        chi2 = 0;
        for (int i = 0; (i < N); i++)
        {
            if (stats->dy[i] > 0)
            {
                dy = stats->dy[i];
            }
            else
            {
                dy = 1;
            }
            chi2 += gmx::square((stats->y[i] - (stats->a * stats->x[i] + stats->b)) / dy);
        }
        if (N > 2)
        {
            stats->chi2 = std::sqrt(chi2 / (N - 2));

            /* Look up equations! */
            dx2            = (xx - sx * sx);
            stats->sigma_a = std::sqrt(stats->chi2 / ((N - 2) * dx2));
            stats->sigma_b = stats->sigma_a * std::sqrt(xx);
            stats->Rfit    = std::abs(ssxy) / std::sqrt(ssxx * ssyy);
        }
        else
        {
            stats->chi2    = 0;
            stats->sigma_a = 0;
            stats->sigma_b = 0;
            stats->Rfit    = 0;
        }

        stats->computed = 1;
    }
}

void gmx_stats_get_ab(gmx_stats_t gstats, int weight, real* a, real* b, real* da, real* db, real* chi2, real* Rfit)
{
    gmx_stats* stats = gstats;

    gmx_stats_compute(stats, weight);
    if (nullptr != a)
    {
        *a = stats->a;
    }
    if (nullptr != b)
    {
        *b = stats->b;
    }
    if (nullptr != da)
    {
        *da = stats->sigma_a;
    }
    if (nullptr != db)
    {
        *db = stats->sigma_b;
    }
    if (nullptr != chi2)
    {
        *chi2 = stats->chi2;
    }
    if (nullptr != Rfit)
    {
        *Rfit = stats->Rfit;
    }
}

real gmx_stats_get_average(gmx_stats_t gstats)
{
    gmx_stats* stats = gstats;

    if (gstats->np < 1)
    {
        GMX_THROW(gmx::InconsistentInputError("No points to average"));
    }
    gmx_stats_compute(stats, elsqWEIGHT_NONE);

    return stats->aver;
}

std::tuple<real, real, real> gmx_stats_get_ase(gmx_stats_t gstats)
{
    gmx_stats* stats = gstats;

    if (gstats->np < 1)
    {
        GMX_THROW(gmx::InconsistentInputError("No points to average"));
    }
    gmx_stats_compute(stats, elsqWEIGHT_NONE);

    return std::make_tuple(stats->aver, stats->sigma_aver, stats->error);
}

// When using GMX_DOUBLE=OFF, some callers want to analyse x values
// that are already computed in double precision. So we need to
// compile two versions, so that the promotion to double is used when
// needed.
template<typename RealT>
void low_lsq_y_ax_b(int n, const RealT* xr, real yr[], real* a, real* b, real* r, real* chi2)
{
    gmx_stats_t lsq = gmx_stats_init();
    for (int i = 0; (i < n); i++)
    {
        gmx_stats_add_point(lsq, double(xr[i]), yr[i], 0, 0);
    }
    gmx_stats_get_ab(lsq, elsqWEIGHT_NONE, a, b, nullptr, nullptr, chi2, r);
    gmx_stats_free(lsq);
}

void lsq_y_ax_b(int n, real x[], real y[], real* a, real* b, real* r, real* chi2)
{
    low_lsq_y_ax_b(n, x, y, a, b, r, chi2);
}

void lsq_y_ax_b_xdouble(int n, double x[], real y[], real* a, real* b, real* r, real* chi2)
{
    low_lsq_y_ax_b(n, x, y, a, b, r, chi2);
}

void lsq_y_ax_b_error(int n, real x[], real y[], real dy[], real* a, real* b, real* da, real* db, real* r, real* chi2)
{
    if (n < 1)
    {
        GMX_THROW(gmx::InconsistentInputError("No points to fit"));
    }

    gmx_stats_t lsq = gmx_stats_init();
    for (int i = 0; (i < n); i++)
    {
        gmx_stats_add_point(lsq, x[i], y[i], 0, dy[i]);
    }
    gmx_stats_get_ab(lsq, elsqWEIGHT_Y, a, b, da, db, chi2, r);
    gmx_stats_free(lsq);
}
