/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include <math.h>
#include "gromacs/legacyheaders/network.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vec.h"
#include "enerdata.h"

typedef struct {
    gmx_int64_t     np;
    double          sum;
    double          sav;
    double          sav2;
} ee_sum_t;

typedef struct {
    int             b;
    ee_sum_t        sum;
    gmx_int64_t     nst;
    gmx_int64_t     nst_min;
} ener_ee_t;

static void clear_ee_sum(ee_sum_t *ees)
{
    ees->sav  = 0;
    ees->sav2 = 0;
    ees->np   = 0;
    ees->sum  = 0;
}

static void add_ee_sum(ee_sum_t *ees, double sum, int np)
{
    ees->np  += np;
    ees->sum += sum;
}

static void add_ee_av(ee_sum_t *ees)
{
    double av;

    av         = ees->sum/ees->np;
    ees->sav  += av;
    ees->sav2 += av*av;
    ees->np    = 0;
    ees->sum   = 0;
}

static double calc_ee2(int nb, ee_sum_t *ees)
{
    return (ees->sav2/nb - dsqr(ees->sav/nb))/(nb - 1);
}

static void set_ee_av(ener_ee_t *eee)
{
    if (debug)
    {
        char buf[STEPSTRSIZE];
        fprintf(debug, "Storing average for err.est.: %s steps\n",
                gmx_step_str(eee->nst, buf));
    }
    add_ee_av(&eee->sum);
    eee->b++;
    if (eee->b == 1 || eee->nst < eee->nst_min)
    {
        eee->nst_min = eee->nst;
    }
    eee->nst = 0;
}

char *ee_pr(double ee, char *buf)
{
    char   tmp[100];
    double rnd;

    if (ee < 0)
    {
        sprintf(buf, "%s", "--");
    }
    else
    {
        /* Round to two decimals by printing. */
        sprintf(tmp, "%12.1e", ee);
        sscanf(tmp, "%12lf", &rnd);
        sprintf(buf, "%g", rnd);
    }

    return buf;
}

void calc_averages(int nset, enerdata_t *edat, int nbmin, int nbmax)
{
    int             nb, i, f, nee;
    double          sum, sum2, sump, see2;
    gmx_int64_t     np, p, bound_nb;
    enerdat_t      *ed;
    exactsum_t     *es;
    gmx_bool        bAllZero;
    double          x, sx, sy, sxx, sxy;
    ener_ee_t      *eee;

    /* Check if we have exact statistics over all points */
    for (i = 0; i < nset; i++)
    {
        ed             = &edat->s[i];
        ed->bExactStat = FALSE;
        if (edat->npoints > 0)
        {
            /* All energy file sum entries 0 signals no exact sums.
             * But if all energy values are 0, we still have exact sums.
             */
            bAllZero = TRUE;
            for (f = 0; f < edat->nframes && !ed->bExactStat; f++)
            {
                if (ed->ener[i] != 0)
                {
                    bAllZero = FALSE;
                }
                ed->bExactStat = (ed->es[f].sum != 0);
            }
            if (bAllZero)
            {
                ed->bExactStat = TRUE;
            }
        }
    }

    snew(eee, nbmax+1);
    for (i = 0; i < nset; i++)
    {
        ed = &edat->s[i];

        sum  = 0;
        sum2 = 0;
        np   = 0;
        sx   = 0;
        sy   = 0;
        sxx  = 0;
        sxy  = 0;
        for (nb = nbmin; nb <= nbmax; nb++)
        {
            eee[nb].b     = 0;
            clear_ee_sum(&eee[nb].sum);
            eee[nb].nst     = 0;
            eee[nb].nst_min = 0;
        }
        for (f = 0; f < edat->nframes; f++)
        {
            es = &ed->es[f];

            if (ed->bExactStat)
            {
                /* Add the sum and the sum of variances to the totals. */
                p     = edat->points[f];
                sump  = es->sum;
                sum2 += es->sum2;
                if (np > 0)
                {
                    sum2 += dsqr(sum/np - (sum + es->sum)/(np + p))
                        *np*(np + p)/p;
                }
            }
            else
            {
                /* Add a single value to the sum and sum of squares. */
                p     = 1;
                sump  = ed->ener[f];
                sum2 += dsqr(sump);
            }

            /* sum has to be increased after sum2 */
            np  += p;
            sum += sump;

            /* For the linear regression use variance 1/p.
             * Note that sump is the sum, not the average, so we don't need p*.
             */
            x    = edat->step[f] - 0.5*(edat->steps[f] - 1);
            sx  += p*x;
            sy  += sump;
            sxx += p*x*x;
            sxy += x*sump;

            for (nb = nbmin; nb <= nbmax; nb++)
            {
                /* Check if the current end step is closer to the desired
                 * block boundary than the next end step.
                 */
                bound_nb = (edat->step[0]-1)*nb + edat->nsteps*(eee[nb].b+1);
                if (eee[nb].nst > 0 &&
                    bound_nb - edat->step[f-1]*nb < edat->step[f]*nb - bound_nb)
                {
                    set_ee_av(&eee[nb]);
                }
                if (f == 0)
                {
                    eee[nb].nst = 1;
                }
                else
                {
                    eee[nb].nst += edat->step[f] - edat->step[f-1];
                }
                if (ed->bExactStat)
                {
                    add_ee_sum(&eee[nb].sum, es->sum, edat->points[f]);
                }
                else
                {
                    add_ee_sum(&eee[nb].sum, edat->s[i].ener[f], 1);
                }
                bound_nb = (edat->step[0]-1)*nb + edat->nsteps*(eee[nb].b+1);
                if (edat->step[f]*nb >= bound_nb)
                {
                    set_ee_av(&eee[nb]);
                }
            }
        }

        edat->s[i].av = sum/np;
        if (ed->bExactStat)
        {
            edat->s[i].rmsd = sqrt(sum2/np);
        }
        else
        {
            edat->s[i].rmsd = sqrt(sum2/np - dsqr(edat->s[i].av));
        }

        if (edat->nframes > 1)
        {
            edat->s[i].slope = (np*sxy - sx*sy)/(np*sxx - sx*sx);
        }
        else
        {
            edat->s[i].slope = 0;
        }

        nee  = 0;
        see2 = 0;
        for (nb = nbmin; nb <= nbmax; nb++)
        {
            /* Check if we actually got nb blocks and if the smallest
             * block is not shorter than 80% of the average.
             */
            if (debug)
            {
                char buf1[STEPSTRSIZE], buf2[STEPSTRSIZE];
                fprintf(debug, "Requested %d blocks, we have %d blocks, min %s nsteps %s\n",
                        nb, eee[nb].b,
                        gmx_step_str(eee[nb].nst_min, buf1),
                        gmx_step_str(edat->nsteps, buf2));
            }
            if (eee[nb].b == nb && 5*nb*eee[nb].nst_min >= 4*edat->nsteps)
            {
                see2 += calc_ee2(nb, &eee[nb].sum);
                nee++;
            }
        }
        if (nee > 0)
        {
            edat->s[i].ee = sqrt(see2/nee);
        }
        else
        {
            edat->s[i].ee = -1;
        }
    }
    sfree(eee);
}

enerdata_t *calc_sum(int nset, enerdata_t *edat, int nbmin, int nbmax)
{
    enerdata_t *esum;
    enerdat_t  *s;
    int         f, i;
    double      sum;

    snew(esum, 1);
    *esum = *edat;
    snew(esum->s, 1);
    s = &esum->s[0];
    snew(s->ener, esum->nframes);
    snew(s->es, esum->nframes);

    s->bExactStat = TRUE;
    s->slope      = 0;
    for (i = 0; i < nset; i++)
    {
        if (!edat->s[i].bExactStat)
        {
            s->bExactStat = FALSE;
        }
        s->slope += edat->s[i].slope;
    }

    for (f = 0; f < edat->nframes; f++)
    {
        sum = 0;
        for (i = 0; i < nset; i++)
        {
            sum += edat->s[i].ener[f];
        }
        s->ener[f] = sum;
        sum        = 0;
        for (i = 0; i < nset; i++)
        {
            sum += edat->s[i].es[f].sum;
        }
        s->es[f].sum  = sum;
        s->es[f].sum2 = 0;
    }

    calc_averages(1, esum, nbmin, nbmax);

    return esum;
}

void remove_drift(int nset, int nbmin, int nbmax, real dt, enerdata_t *edat)
{
/* Remove the drift by performing a fit to y = ax+b.
   Uses 5 iterations. */
    int    i, j, k;
    double delta;

    edat->npoints = edat->nframes;
    edat->nsteps  = edat->nframes;

    for (k = 0; (k < 5); k++)
    {
        for (i = 0; (i < nset); i++)
        {
            delta = edat->s[i].slope*dt;

            if (NULL != debug)
            {
                fprintf(debug, "slope for set %d is %g\n", i, edat->s[i].slope);
            }

            for (j = 0; (j < edat->nframes); j++)
            {
                edat->s[i].ener[j]   -= j*delta;
                edat->s[i].es[j].sum  = 0;
                edat->s[i].es[j].sum2 = 0;
            }
        }
        calc_averages(nset, edat, nbmin, nbmax);
    }
}
