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
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"


namespace gmx
{

struct stats
{
    double  aa, a, b, sigma_a, sigma_b, aver, sigma_aver, error;
    double  rmsd, Rfit, chi2;
    double *x, *y, *dx, *dy;
    int     computed;
    int     np, nalloc;
};

stats_t stats_init()
{
    stats* s;

    snew(s, 1);

    return static_cast<stats_t>(s);
}

void stats_free(stats_t gstats)
{
    stats* s = static_cast<stats*>(gstats);

    sfree(s->x);
    sfree(s->y);
    sfree(s->dx);
    sfree(s->dy);
    sfree(s);
}

void stats_add_point(stats_t gstats, double x, double y, double dx, double dy)
{
    stats* s = gstats;

    if (s->np + 1 >= s->nalloc)
    {
        if (s->nalloc == 0)
        {
            s->nalloc = 1024;
        }
        else
        {
            s->nalloc *= 2;
        }
        srenew(s->x, s->nalloc);
        srenew(s->y, s->nalloc);
        srenew(s->dx, s->nalloc);
        srenew(s->dy, s->nalloc);
        for (int i = s->np; (i < s->nalloc); i++)
        {
            s->x[i]  = 0;
            s->y[i]  = 0;
            s->dx[i] = 0;
            s->dy[i] = 0;
        }
    }
    s->x[s->np]  = x;
    s->y[s->np]  = y;
    s->dx[s->np] = dx;
    s->dy[s->np] = dy;
    s->np++;
    s->computed = 0;
}

static void stats_compute(stats* s, int weight)
{
    double yx, xx, sx, sy, dy, chi2, d2;
    double ssxx, ssyy, ssxy;
    double w, wtot, yx_nw, sy_nw, sx_nw, yy_nw, xx_nw, dx2;

    int N = s->np;

    if (s->computed == 0)
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
            d2 += square(s->x[i] - s->y[i]);
            if (((s->dy[i]) != 0.0) && (weight == elsqWEIGHT_Y))
            {
                w = 1 / square(s->dy[i]);
            }
            else
            {
                w = 1;
            }

            wtot += w;

            xx += w * square(s->x[i]);
            xx_nw += square(s->x[i]);

            yy_nw += square(s->y[i]);

            yx += w * s->y[i] * s->x[i];
            yx_nw += s->y[i] * s->x[i];

            sx += w * s->x[i];
            sx_nw += s->x[i];

            sy += w * s->y[i];
            sy_nw += s->y[i];
        }

        /* Compute average, sigma and error */
        s->aver       = sy_nw / N;
        s->sigma_aver = std::sqrt(yy_nw / N - square(sy_nw / N));
        s->error      = s->sigma_aver / std::sqrt(static_cast<double>(N));

        /* Compute RMSD between x and y */
        s->rmsd = std::sqrt(d2 / N);

        /* Correlation coefficient for data */
        yx_nw /= N;
        xx_nw /= N;
        yy_nw /= N;
        sx_nw /= N;
        sy_nw /= N;
        ssxx = N * (xx_nw - square(sx_nw));
        ssyy = N * (yy_nw - square(sy_nw));
        ssxy = N * (yx_nw - (sx_nw * sy_nw));

        /* Compute straight line through datapoints, either with intercept
           zero (result in aa) or with intercept variable (results in a
           and b) */
        yx = yx / wtot;
        xx = xx / wtot;
        sx = sx / wtot;
        sy = sy / wtot;

        s->aa = (yx / xx);
        s->a  = (yx - sx * sy) / (xx - sx * sx);
        s->b  = (sy) - (s->a) * (sx);

        /* Compute chi2, deviation from a line y = ax+b. */
        chi2 = 0;
        for (int i = 0; (i < N); i++)
        {
            if (s->dy[i] > 0)
            {
                dy = s->dy[i];
            }
            else
            {
                dy = 1;
            }
            chi2 += square((s->y[i] - (s->a * s->x[i] + s->b)) / dy);
        }
        if (N > 2)
        {
            s->chi2 = std::sqrt(chi2 / (N - 2));

            /* Look up equations! */
            dx2        = (xx - sx * sx);
            s->sigma_a = std::sqrt(s->chi2 / ((N - 2) * dx2));
            s->sigma_b = s->sigma_a * std::sqrt(xx);
            s->Rfit    = std::abs(ssxy) / std::sqrt(ssxx * ssyy);
        }
        else
        {
            s->chi2    = 0;
            s->sigma_a = 0;
            s->sigma_b = 0;
            s->Rfit    = 0;
        }

        s->computed = 1;
    }
}

void stats_get_ab(stats_t gstats, int weight, real* a, real* b, real* da, real* db, real* chi2, real* Rfit)
{
    stats* s = gstats;

    stats_compute(s, weight);
    if (nullptr != a)
    {
        *a = s->a;
    }
    if (nullptr != b)
    {
        *b = s->b;
    }
    if (nullptr != da)
    {
        *da = s->sigma_a;
    }
    if (nullptr != db)
    {
        *db = s->sigma_b;
    }
    if (nullptr != chi2)
    {
        *chi2 = s->chi2;
    }
    if (nullptr != Rfit)
    {
        *Rfit = s->Rfit;
    }
}

real stats_get_average(stats_t gstats)
{
    stats* s = gstats;

    if (s->np < 1)
    {
        GMX_THROW(InconsistentInputError("No points to average"));
    }
    stats_compute(s, elsqWEIGHT_NONE);

    return s->aver;
}

std::tuple<real, real, real> stats_get_ase(stats_t gstats)
{
    stats* s = gstats;

    if (s->np < 1)
    {
        GMX_THROW(InconsistentInputError("No points to average"));
    }
    stats_compute(s, elsqWEIGHT_NONE);

    return std::make_tuple(s->aver, s->sigma_aver, s->error);
}

// When using GMX_DOUBLE=OFF, some callers want to analyse x values
// that are already computed in double precision. So we need to
// compile two versions, so that the promotion to double is used when
// needed.
template<typename RealT>
void low_lsq_y_ax_b(int n, const RealT* xr, real yr[], real* a, real* b, real* r, real* chi2)
{
    stats_t lsq = stats_init();
    for (int i = 0; (i < n); i++)
    {
        stats_add_point(lsq, double(xr[i]), yr[i], 0, 0);
    }
    stats_get_ab(lsq, elsqWEIGHT_NONE, a, b, nullptr, nullptr, chi2, r);
    stats_free(lsq);
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
        GMX_THROW(InconsistentInputError("No points to fit"));
    }

    stats_t lsq = stats_init();
    for (int i = 0; (i < n); i++)
    {
        stats_add_point(lsq, x[i], y[i], 0, dy[i]);
    }
    stats_get_ab(lsq, elsqWEIGHT_Y, a, b, da, db, chi2, r);
    stats_free(lsq);
}

} // namespace gmx
