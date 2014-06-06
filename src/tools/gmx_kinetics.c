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
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "gmx_fatal.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "readinp.h"
#include "statutil.h"
#include "txtdump.h"
#include "gstat.h"
#include "xvgr.h"
#include "physics.h"
#include "gmx_ana.h"

enum {
    epAuf, epEuf, epAfu, epEfu, epNR
};
enum {
    eqAif, eqEif, eqAfi, eqEfi, eqAui, eqEui, eqAiu, eqEiu, eqNR
};
static char *eep[epNR] = { "Af", "Ef", "Au", "Eu" };
static char *eeq[eqNR] = { "Aif", "Eif", "Afi", "Efi", "Aui", "Eui", "Aiu", "Eiu" };

typedef struct {
    int       nreplica;  /* Number of replicas in the calculation                   */
    int       nframe;    /* Number of time frames                                   */
    int       nstate;    /* Number of states the system can be in, e.g. F,I,U       */
    int       nparams;   /* Is 2, 4 or 8                                            */
    gmx_bool *bMask;     /* Determine whether this replica is part of the d2 comp.  */
    gmx_bool  bSum;
    gmx_bool  bDiscrete; /* Use either discrete folding (0/1) or a continuous       */
    /* criterion */
    int       nmask;     /* Number of replicas taken into account                   */
    real      dt;        /* Timestep between frames                                 */
    int       j0, j1;    /* Range of frames used in calculating delta               */
    real    **temp, **data, **data2;
    int     **state;     /* State index running from 0 (F) to nstate-1 (U)          */
    real    **beta, **fcalt, **icalt;
    real     *time, *sumft, *sumit, *sumfct, *sumict;
    real     *params;
    real     *d2_replica;
} t_remd_data;

#ifdef HAVE_LIBGSL
#include <gsl/gsl_multimin.h>

static char *itoa(int i)
{
    static char ptr[12];

    sprintf(ptr, "%d", i);
    return ptr;
}

static char *epnm(int nparams, int index)
{
    static char buf[32], from[8], to[8];
    int         nn, ni, ii;

    range_check(index, 0, nparams);
    if ((nparams == 2) || (nparams == 4))
    {
        return eep[index];
    }
    else if ((nparams > 4) && (nparams % 4 == 0))
    {
        return eeq[index];
    }
    else
    {
        gmx_fatal(FARGS, "Don't know how to handle %d parameters", nparams);
    }

    return NULL;
}

static gmx_bool bBack(t_remd_data *d)
{
    return (d->nparams > 2);
}

static real is_folded(t_remd_data *d, int irep, int jframe)
{
    if (d->state[irep][jframe] == 0)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

static real is_unfolded(t_remd_data *d, int irep, int jframe)
{
    if (d->state[irep][jframe] == d->nstate-1)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

static real is_intermediate(t_remd_data *d, int irep, int jframe)
{
    if ((d->state[irep][jframe] == 1) && (d->nstate > 2))
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

static void integrate_dfdt(t_remd_data *d)
{
    int    i, j;
    double beta, ddf, ddi, df, db, fac, sumf, sumi, area;

    d->sumfct[0] = 0;
    d->sumict[0] = 0;
    for (i = 0; (i < d->nreplica); i++)
    {
        if (d->bMask[i])
        {
            if (d->bDiscrete)
            {
                ddf = 0.5*d->dt*is_folded(d, i, 0);
                ddi = 0.5*d->dt*is_intermediate(d, i, 0);
            }
            else
            {
                ddf = 0.5*d->dt*d->state[i][0];
                ddi = 0.0;
            }
            d->fcalt[i][0] = ddf;
            d->icalt[i][0] = ddi;
            d->sumfct[0]  += ddf;
            d->sumict[0]  += ddi;
        }
    }
    for (j = 1; (j < d->nframe); j++)
    {
        if (j == d->nframe-1)
        {
            fac = 0.5*d->dt;
        }
        else
        {
            fac = d->dt;
        }
        sumf = sumi = 0;
        for (i = 0; (i < d->nreplica); i++)
        {
            if (d->bMask[i])
            {
                beta = d->beta[i][j];
                if ((d->nstate <= 2) || d->bDiscrete)
                {
                    if (d->bDiscrete)
                    {
                        df = (d->params[epAuf]*exp(-beta*d->params[epEuf])*
                              is_unfolded(d, i, j));
                    }
                    else
                    {
                        area = (d->data2 ? d->data2[i][j] : 1.0);
                        df   =  area*d->params[epAuf]*exp(-beta*d->params[epEuf]);
                    }
                    if (bBack(d))
                    {
                        db = 0;
                        if (d->bDiscrete)
                        {
                            db = (d->params[epAfu]*exp(-beta*d->params[epEfu])*
                                  is_folded(d, i, j));
                        }
                        else
                        {
                            gmx_fatal(FARGS, "Back reaction not implemented with continuous");
                        }
                        ddf = fac*(df-db);
                    }
                    else
                    {
                        ddf = fac*df;
                    }
                    d->fcalt[i][j] = d->fcalt[i][j-1] + ddf;
                    sumf          += ddf;
                }
                else
                {
                    ddf = fac*((d->params[eqAif]*exp(-beta*d->params[eqEif])*
                                is_intermediate(d, i, j)) -
                               (d->params[eqAfi]*exp(-beta*d->params[eqEfi])*
                                is_folded(d, i, j)));
                    ddi = fac*((d->params[eqAui]*exp(-beta*d->params[eqEui])*
                                is_unfolded(d, i, j)) -
                               (d->params[eqAiu]*exp(-beta*d->params[eqEiu])*
                                is_intermediate(d, i, j)));
                    d->fcalt[i][j] = d->fcalt[i][j-1] + ddf;
                    d->icalt[i][j] = d->icalt[i][j-1] + ddi;
                    sumf          += ddf;
                    sumi          += ddi;
                }
            }
        }
        d->sumfct[j] = d->sumfct[j-1] + sumf;
        d->sumict[j] = d->sumict[j-1] + sumi;
    }
    if (debug)
    {
        fprintf(debug, "@type xy\n");
        for (j = 0; (j < d->nframe); j++)
        {
            fprintf(debug, "%8.3f  %12.5e\n", d->time[j], d->sumfct[j]);
        }
        fprintf(debug, "&\n");
    }
}

static void sum_ft(t_remd_data *d)
{
    int    i, j;
    double fac;

    for (j = 0; (j < d->nframe); j++)
    {
        d->sumft[j] = 0;
        d->sumit[j] = 0;
        if ((j == 0) || (j == d->nframe-1))
        {
            fac = d->dt*0.5;
        }
        else
        {
            fac = d->dt;
        }
        for (i = 0; (i < d->nreplica); i++)
        {
            if (d->bMask[i])
            {
                if (d->bDiscrete)
                {
                    d->sumft[j] += fac*is_folded(d, i, j);
                    d->sumit[j] += fac*is_intermediate(d, i, j);
                }
                else
                {
                    d->sumft[j] += fac*d->state[i][j];
                }
            }
        }
    }
}

static double calc_d2(t_remd_data *d)
{
    int    i, j;
    double dd2, d2 = 0, dr2, tmp;

    integrate_dfdt(d);

    if (d->bSum)
    {
        for (j = d->j0; (j < d->j1); j++)
        {
            if (d->bDiscrete)
            {
                d2  += sqr(d->sumft[j]-d->sumfct[j]);
                if (d->nstate > 2)
                {
                    d2 += sqr(d->sumit[j]-d->sumict[j]);
                }
            }
            else
            {
                d2  += sqr(d->sumft[j]-d->sumfct[j]);
            }
        }
    }
    else
    {
        for (i = 0; (i < d->nreplica); i++)
        {
            dr2 = 0;
            if (d->bMask[i])
            {
                for (j = d->j0; (j < d->j1); j++)
                {
                    tmp  = sqr(is_folded(d, i, j)-d->fcalt[i][j]);
                    d2  += tmp;
                    dr2 += tmp;
                    if (d->nstate > 2)
                    {
                        tmp  = sqr(is_intermediate(d, i, j)-d->icalt[i][j]);
                        d2  += tmp;
                        dr2 += tmp;
                    }
                }
                d->d2_replica[i] = dr2/(d->j1-d->j0);
            }
        }
    }
    dd2 = (d2/(d->j1-d->j0))/(d->bDiscrete ? d->nmask : 1);

    return dd2;
}

static double my_f(const gsl_vector *v, void *params)
{
    t_remd_data *d       = (t_remd_data *) params;
    double       penalty = 0;
    int          i;

    for (i = 0; (i < d->nparams); i++)
    {
        d->params[i] = gsl_vector_get(v, i);
        if (d->params[i] < 0)
        {
            penalty += 10;
        }
    }
    if (penalty > 0)
    {
        return penalty;
    }
    else
    {
        return calc_d2(d);
    }
}

static void optimize_remd_parameters(FILE *fp, t_remd_data *d, int maxiter,
                                     real tol)
{
    real   size, d2;
    int    iter   = 0;
    int    status = 0;
    int    i;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer            *s;

    gsl_vector                         *x, *dx;
    gsl_multimin_function               my_func;

    my_func.f      = &my_f;
    my_func.n      = d->nparams;
    my_func.params = (void *) d;

    /* Starting point */
    x = gsl_vector_alloc (my_func.n);
    for (i = 0; (i < my_func.n); i++)
    {
        gsl_vector_set (x, i, d->params[i]);
    }

    /* Step size, different for each of the parameters */
    dx = gsl_vector_alloc (my_func.n);
    for (i = 0; (i < my_func.n); i++)
    {
        gsl_vector_set (dx, i, 0.1*d->params[i]);
    }

    T = gsl_multimin_fminimizer_nmsimplex;
    s = gsl_multimin_fminimizer_alloc (T, my_func.n);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);
    gsl_vector_free (x);
    gsl_vector_free (dx);

    printf ("%5s", "Iter");
    for (i = 0; (i < my_func.n); i++)
    {
        printf(" %12s", epnm(my_func.n, i));
    }
    printf (" %12s %12s\n", "NM Size", "Chi2");

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);

        if (status != 0)
        {
            gmx_fatal(FARGS, "Something went wrong in the iteration in minimizer %s",
                      gsl_multimin_fminimizer_name(s));
        }

        d2     = gsl_multimin_fminimizer_minimum(s);
        size   = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, tol);

        if (status == GSL_SUCCESS)
        {
            printf ("Minimum found using %s at:\n",
                    gsl_multimin_fminimizer_name(s));
        }

        printf ("%5d", iter);
        for (i = 0; (i < my_func.n); i++)
        {
            printf(" %12.4e", gsl_vector_get (s->x, i));
        }
        printf (" %12.4e %12.4e\n", size, d2);
    }
    while ((status == GSL_CONTINUE) && (iter < maxiter));

    gsl_multimin_fminimizer_free (s);
}

static void preprocess_remd(FILE *fp, t_remd_data *d, real cutoff, real tref,
                            real ucut, gmx_bool bBack, real Euf, real Efu,
                            real Ei, real t0, real t1, gmx_bool bSum, gmx_bool bDiscrete,
                            int nmult)
{
    int  i, j, ninter;
    real dd, tau_f, tau_u;

    ninter = (ucut > cutoff) ? 1 : 0;
    if (ninter && (ucut <= cutoff))
    {
        gmx_fatal(FARGS, "You have requested an intermediate but the cutoff for intermediates %f is smaller than the normal cutoff(%f)", ucut, cutoff);
    }

    if (!bBack)
    {
        d->nparams = 2;
        d->nstate  = 2;
    }
    else
    {
        d->nparams = 4*(1+ninter);
        d->nstate  = 2+ninter;
    }
    d->bSum      = bSum;
    d->bDiscrete = bDiscrete;
    snew(d->beta, d->nreplica);
    snew(d->state, d->nreplica);
    snew(d->bMask, d->nreplica);
    snew(d->d2_replica, d->nreplica);
    snew(d->sumft, d->nframe);
    snew(d->sumit, d->nframe);
    snew(d->sumfct, d->nframe);
    snew(d->sumict, d->nframe);
    snew(d->params, d->nparams);
    snew(d->fcalt, d->nreplica);
    snew(d->icalt, d->nreplica);

    /* convert_times(d->nframe,d->time); */

    if (t0 < 0)
    {
        d->j0 = 0;
    }
    else
    {
        for (d->j0 = 0; (d->j0 < d->nframe) && (d->time[d->j0] < t0); d->j0++)
        {
            ;
        }
    }
    if (t1 < 0)
    {
        d->j1 = d->nframe;
    }
    else
    {
        for (d->j1 = 0; (d->j1 < d->nframe) && (d->time[d->j1] < t1); d->j1++)
        {
            ;
        }
    }
    if ((d->j1-d->j0) < d->nparams+2)
    {
        gmx_fatal(FARGS, "Start (%f) or end time (%f) for fitting inconsistent. Reduce t0, increase t1 or supply more data", t0, t1);
    }
    fprintf(fp, "Will optimize from %g to %g\n",
            d->time[d->j0], d->time[d->j1-1]);
    d->nmask = d->nreplica;
    for (i = 0; (i < d->nreplica); i++)
    {
        snew(d->beta[i], d->nframe);
        snew(d->state[i], d->nframe);
        snew(d->fcalt[i], d->nframe);
        snew(d->icalt[i], d->nframe);
        d->bMask[i] = TRUE;
        for (j = 0; (j < d->nframe); j++)
        {
            d->beta[i][j] = 1.0/(BOLTZ*d->temp[i][j]);
            dd            = d->data[i][j];
            if (bDiscrete)
            {
                if (dd <= cutoff)
                {
                    d->state[i][j] = 0;
                }
                else if ((ucut > cutoff) && (dd <= ucut))
                {
                    d->state[i][j] = 1;
                }
                else
                {
                    d->state[i][j] = d->nstate-1;
                }
            }
            else
            {
                d->state[i][j] = dd*nmult;
            }
        }
    }
    sum_ft(d);

    /* Assume forward rate constant is half the total time in this
     * simulation and backward is ten times as long */
    if (bDiscrete)
    {
        tau_f            = d->time[d->nframe-1];
        tau_u            = 4*tau_f;
        d->params[epEuf] = Euf;
        d->params[epAuf] = exp(d->params[epEuf]/(BOLTZ*tref))/tau_f;
        if (bBack)
        {
            d->params[epEfu] = Efu;
            d->params[epAfu] = exp(d->params[epEfu]/(BOLTZ*tref))/tau_u;
            if (ninter > 0)
            {
                d->params[eqEui] = Ei;
                d->params[eqAui] = exp(d->params[eqEui]/(BOLTZ*tref))/tau_u;
                d->params[eqEiu] = Ei;
                d->params[eqAiu] = exp(d->params[eqEiu]/(BOLTZ*tref))/tau_u;
            }
        }
        else
        {
            d->params[epAfu]  = 0;
            d->params[epEfu]  = 0;
        }
    }
    else
    {
        d->params[epEuf] = Euf;
        if (d->data2)
        {
            d->params[epAuf] = 0.1;
        }
        else
        {
            d->params[epAuf] = 20.0;
        }
    }
}

static real tau(real A, real E, real T)
{
    return exp(E/(BOLTZ*T))/A;
}

static real folded_fraction(t_remd_data *d, real tref)
{
    real tauf, taub;

    tauf = tau(d->params[epAuf], d->params[epEuf], tref);
    taub = tau(d->params[epAfu], d->params[epEfu], tref);

    return (taub/(tauf+taub));
}

static void print_tau(FILE *gp, t_remd_data *d, real tref)
{
    real tauf, taub, ddd, fff, DG, DH, TDS, Tm, Tb, Te, Fb, Fe, Fm;
    int  i, np = d->nparams;

    ddd = calc_d2(d);
    fprintf(gp, "Final value for Chi2 = %12.5e (%d replicas)\n", ddd, d->nmask);
    tauf = tau(d->params[epAuf], d->params[epEuf], tref);
    fprintf(gp, "%s = %12.5e %s = %12.5e (kJ/mole)\n",
            epnm(np, epAuf), d->params[epAuf],
            epnm(np, epEuf), d->params[epEuf]);
    if (bBack(d))
    {
        taub = tau(d->params[epAfu], d->params[epEfu], tref);
        fprintf(gp, "%s = %12.5e %s = %12.5e (kJ/mole)\n",
                epnm(np, epAfu), d->params[epAfu],
                epnm(np, epEfu), d->params[epEfu]);
        fprintf(gp, "Equilibrium properties at T = %g\n", tref);
        fprintf(gp, "tau_f = %8.3f ns, tau_b = %8.3f ns\n", tauf/1000, taub/1000);
        fff = taub/(tauf+taub);
        DG  = BOLTZ*tref*log(fff/(1-fff));
        DH  = d->params[epEfu]-d->params[epEuf];
        TDS = DH-DG;
        fprintf(gp, "Folded fraction     F = %8.3f\n", fff);
        fprintf(gp, "Unfolding energies: DG = %8.3f  DH = %8.3f TDS = %8.3f\n",
                DG, DH, TDS);
        Tb = 260;
        Te = 420;
        Tm = 0;
        Fm = 0;
        Fb = folded_fraction(d, Tb);
        Fe = folded_fraction(d, Te);
        while ((Te-Tb > 0.001) && (Fm != 0.5))
        {
            Tm = 0.5*(Tb+Te);
            Fm = folded_fraction(d, Tm);
            if (Fm > 0.5)
            {
                Fb = Fm;
                Tb = Tm;
            }
            else if (Fm < 0.5)
            {
                Te = Tm;
                Fe = Fm;
            }
        }
        if ((Fb-0.5)*(Fe-0.5) <= 0)
        {
            fprintf(gp, "Melting temperature Tm = %8.3f K\n", Tm);
        }
        else
        {
            fprintf(gp, "No melting temperature detected between 260 and 420K\n");
        }
        if (np > 4)
        {
            char *ptr;
            fprintf(gp, "Data for intermediates at T = %g\n", tref);
            fprintf(gp, "%8s  %10s  %10s  %10s\n", "Name", "A", "E", "tau");
            for (i = 0; (i < np/2); i++)
            {
                tauf = tau(d->params[2*i], d->params[2*i+1], tref);
                ptr  = epnm(d->nparams, 2*i);
                fprintf(gp, "%8s  %10.3e  %10.3e  %10.3e\n", ptr+1,
                        d->params[2*i], d->params[2*i+1], tauf/1000);
            }
        }
    }
    else
    {
        fprintf(gp, "Equilibrium properties at T = %g\n", tref);
        fprintf(gp, "tau_f = %8.3f\n", tauf);
    }
}

static void dump_remd_parameters(FILE *gp, t_remd_data *d, const char *fn,
                                 const char *fn2, const char *rfn,
                                 const char *efn, const char *mfn, int skip, real tref,
                                 output_env_t oenv)
{
    FILE       *fp, *hp;
    int         i, j, np = d->nparams;
    real        rhs, tauf, taub, fff, DG;
    real       *params;
    const char *leg[]  = { "Measured", "Fit", "Difference" };
    const char *mleg[] = { "Folded fraction", "DG (kJ/mole)"};
    char      **rleg;
    real        fac[] = { 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03 };
#define NFAC asize(fac)
    real        d2[NFAC];
    double      norm;

    integrate_dfdt(d);
    print_tau(gp, d, tref);
    norm = (d->bDiscrete ? 1.0/d->nmask : 1.0);

    if (fn)
    {
        fp = xvgropen(fn, "Optimized fit to data", "Time (ps)", "Fraction Folded", oenv);
        xvgr_legend(fp, asize(leg), leg, oenv);
        for (i = 0; (i < d->nframe); i++)
        {
            if ((skip <= 0) || ((i % skip) == 0))
            {
                fprintf(fp, "%12.5e  %12.5e  %12.5e  %12.5e\n", d->time[i],
                        d->sumft[i]*norm, d->sumfct[i]*norm,
                        (d->sumft[i]-d->sumfct[i])*norm);
            }
        }
        ffclose(fp);
    }
    if (!d->bSum && rfn)
    {
        snew(rleg, d->nreplica*2);
        for (i = 0; (i < d->nreplica); i++)
        {
            snew(rleg[2*i], 32);
            snew(rleg[2*i+1], 32);
            sprintf(rleg[2*i], "\\f{4}F(t) %d", i);
            sprintf(rleg[2*i+1], "\\f{12}F \\f{4}(t) %d", i);
        }
        fp = xvgropen(rfn, "Optimized fit to data", "Time (ps)", "Fraction Folded", oenv);
        xvgr_legend(fp, d->nreplica*2, (const char**)rleg, oenv);
        for (j = 0; (j < d->nframe); j++)
        {
            if ((skip <= 0) || ((j % skip) == 0))
            {
                fprintf(fp, "%12.5e", d->time[j]);
                for (i = 0; (i < d->nreplica); i++)
                {
                    fprintf(fp, "  %5f  %9.2e", is_folded(d, i, j), d->fcalt[i][j]);
                }
                fprintf(fp, "\n");
            }
        }
        ffclose(fp);
    }

    if (fn2 && (d->nstate > 2))
    {
        fp = xvgropen(fn2, "Optimized fit to data", "Time (ps)",
                      "Fraction Intermediate", oenv);
        xvgr_legend(fp, asize(leg), leg, oenv);
        for (i = 0; (i < d->nframe); i++)
        {
            if ((skip <= 0) || ((i % skip) == 0))
            {
                fprintf(fp, "%12.5e  %12.5e  %12.5e  %12.5e\n", d->time[i],
                        d->sumit[i]*norm, d->sumict[i]*norm,
                        (d->sumit[i]-d->sumict[i])*norm);
            }
        }
        ffclose(fp);
    }
    if (mfn)
    {
        if (bBack(d))
        {
            fp = xvgropen(mfn, "Melting curve", "T (K)", "", oenv);
            xvgr_legend(fp, asize(mleg), mleg, oenv);
            for (i = 260; (i <= 420); i++)
            {
                tauf = tau(d->params[epAuf], d->params[epEuf], 1.0*i);
                taub = tau(d->params[epAfu], d->params[epEfu], 1.0*i);
                fff  = taub/(tauf+taub);
                DG   = BOLTZ*i*log(fff/(1-fff));
                fprintf(fp, "%5d  %8.3f  %8.3f\n", i, fff, DG);
            }
            ffclose(fp);
        }
    }

    if (efn)
    {
        snew(params, d->nparams);
        for (i = 0; (i < d->nparams); i++)
        {
            params[i] = d->params[i];
        }

        hp = xvgropen(efn, "Chi2 as a function of relative parameter",
                      "Fraction", "Chi2", oenv);
        for (j = 0; (j < d->nparams); j++)
        {
            /* Reset all parameters to optimized values */
            if(output_env_get_print_xvgr_codes(oenv))
            {
                fprintf(hp, "@type xy\n");
            }
            for (i = 0; (i < d->nparams); i++)
            {
                d->params[i] = params[i];
            }
            /* Now modify one of them */
            for (i = 0; (i < NFAC); i++)
            {
                d->params[j] = fac[i]*params[j];
                d2[i]        = calc_d2(d);
                fprintf(gp, "%s = %12g  d2 = %12g\n", epnm(np, j), d->params[j], d2[i]);
                fprintf(hp, "%12g  %12g\n", fac[i], d2[i]);
            }
            fprintf(hp, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
        }
        ffclose(hp);
        for (i = 0; (i < d->nparams); i++)
        {
            d->params[i] = params[i];
        }
        sfree(params);
    }
    if (!d->bSum)
    {
        for (i = 0; (i < d->nreplica); i++)
        {
            fprintf(gp, "Chi2[%3d] = %8.2e\n", i, d->d2_replica[i]);
        }
    }
}
#endif /*HAVE_LIBGSL*/

int gmx_kinetics(int argc, char *argv[])
{
    const char     *desc[] = {
        "[TT]g_kinetics[tt] reads two [TT].xvg[tt] files, each one containing data for N replicas.",
        "The first file contains the temperature of each replica at each timestep,",
        "and the second contains real values that can be interpreted as",
        "an indicator for folding. If the value in the file is larger than",
        "the cutoff it is taken to be unfolded and the other way around.[PAR]",
        "From these data an estimate of the forward and backward rate constants",
        "for folding is made at a reference temperature. In addition,",
        "a theoretical melting curve and free energy as a function of temperature",
        "are printed in an [TT].xvg[tt] file.[PAR]",
        "The user can give a max value to be regarded as intermediate",
        "([TT]-ucut[tt]), which, when given will trigger the use of an intermediate state",
        "in the algorithm to be defined as those structures that have",
        "cutoff < DATA < ucut. Structures with DATA values larger than ucut will",
        "not be regarded as potential folders. In this case 8 parameters are optimized.[PAR]",
        "The average fraction foled is printed in an [TT].xvg[tt] file together with the fit to it.",
        "If an intermediate is used a further file will show the build of the intermediate and the fit to that process.[PAR]",
        "The program can also be used with continuous variables (by setting",
        "[TT]-nodiscrete[tt]). In this case kinetics of other processes can be",
        "studied. This is very much a work in progress and hence the manual",
        "(this information) is lagging behind somewhat.[PAR]",
        "In order to compile this program you need access to the GNU",
        "scientific library."
    };
    static int      nreplica  = 1;
    static real     tref      = 298.15;
    static real     cutoff    = 0.2;
    static real     ucut      = 0.0;
    static real     Euf       = 10;
    static real     Efu       = 30;
    static real     Ei        = 10;
    static gmx_bool bHaveT    = TRUE;
    static real     t0        = -1;
    static real     t1        = -1;
    static real     tb        = 0;
    static real     te        = 0;
    static real     tol       = 1e-3;
    static int      maxiter   = 100;
    static int      skip      = 0;
    static int      nmult     = 1;
    static gmx_bool bBack     = TRUE;
    static gmx_bool bSplit    = TRUE;
    static gmx_bool bSum      = TRUE;
    static gmx_bool bDiscrete = TRUE;
    t_pargs         pa[]      = {
        { "-time",    FALSE, etBOOL, {&bHaveT},
          "Expect a time in the input" },
        { "-b",       FALSE, etREAL, {&tb},
          "First time to read from set" },
        { "-e",       FALSE, etREAL, {&te},
          "Last time to read from set" },
        { "-bfit",    FALSE, etREAL, {&t0},
          "Time to start the fit from" },
        { "-efit",    FALSE, etREAL, {&t1},
          "Time to end the fit" },
        { "-T",       FALSE, etREAL, {&tref},
          "Reference temperature for computing rate constants" },
        { "-n",       FALSE, etINT, {&nreplica},
          "Read data for this number of replicas. Only necessary when files are written in xmgrace format using @type and & as delimiters." },
        { "-cut",     FALSE, etREAL, {&cutoff},
          "Cut-off (max) value for regarding a structure as folded" },
        { "-ucut",    FALSE, etREAL, {&ucut},
          "Cut-off (max) value for regarding a structure as intermediate (if not folded)" },
        { "-euf",     FALSE, etREAL, {&Euf},
          "Initial guess for energy of activation for folding (kJ/mol)" },
        { "-efu",     FALSE, etREAL, {&Efu},
          "Initial guess for energy of activation for unfolding (kJ/mol)" },
        { "-ei",      FALSE, etREAL, {&Ei},
          "Initial guess for energy of activation for intermediates (kJ/mol)" },
        { "-maxiter", FALSE, etINT, {&maxiter},
          "Max number of iterations" },
        { "-back",    FALSE, etBOOL, {&bBack},
          "Take the back reaction into account" },
        { "-tol",     FALSE, etREAL, {&tol},
          "Absolute tolerance for convergence of the Nelder and Mead simplex algorithm" },
        { "-skip",    FALSE, etINT, {&skip},
          "Skip points in the output [TT].xvg[tt] file" },
        { "-split",   FALSE, etBOOL, {&bSplit},
          "Estimate error by splitting the number of replicas in two and refitting" },
        { "-sum",     FALSE, etBOOL, {&bSum},
          "Average folding before computing [GRK]chi[grk]^2" },
        { "-discrete", FALSE, etBOOL, {&bDiscrete},
          "Use a discrete folding criterion (F <-> U) or a continuous one" },
        { "-mult",    FALSE, etINT, {&nmult},
          "Factor to multiply the data with before discretization" }
    };
#define NPA asize(pa)

    FILE        *fp;
    real         dt_t, dt_d, dt_d2;
    int          nset_t, nset_d, nset_d2, n_t, n_d, n_d2, i;
    const char  *tfile, *dfile, *dfile2;
    t_remd_data  remd;
    output_env_t oenv;

    t_filenm     fnm[] = {
        { efXVG, "-f",    "temp",    ffREAD   },
        { efXVG, "-d",    "data",    ffREAD   },
        { efXVG, "-d2",   "data2",   ffOPTRD  },
        { efXVG, "-o",    "ft_all",  ffWRITE  },
        { efXVG, "-o2",   "it_all",  ffOPTWR  },
        { efXVG, "-o3",   "ft_repl", ffOPTWR  },
        { efXVG, "-ee",   "err_est", ffOPTWR  },
        { efLOG, "-g",    "remd",    ffWRITE  },
        { efXVG, "-m",    "melt",    ffWRITE  }
    };
#define NFILE asize(fnm)

    CopyRight(stderr, argv[0]);
    parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_BE_NICE | PCA_TIME_UNIT,
                      NFILE, fnm, NPA, pa, asize(desc), desc, 0, NULL, &oenv);

#ifdef HAVE_LIBGSL
    please_cite(stdout, "Spoel2006d");
    if (cutoff < 0)
    {
        gmx_fatal(FARGS, "cutoff should be >= 0 (rather than %f)", cutoff);
    }

    tfile   = opt2fn("-f", NFILE, fnm);
    dfile   = opt2fn("-d", NFILE, fnm);
    dfile2  = opt2fn_null("-d2", NFILE, fnm);

    fp = ffopen(opt2fn("-g", NFILE, fnm), "w");

    remd.temp = read_xvg_time(tfile, bHaveT,
                              opt2parg_bSet("-b", NPA, pa), tb,
                              opt2parg_bSet("-e", NPA, pa), te,
                              nreplica, &nset_t, &n_t, &dt_t, &remd.time);
    printf("Read %d sets of %d points in %s, dt = %g\n\n", nset_t, n_t, tfile, dt_t);
    sfree(remd.time);

    remd.data = read_xvg_time(dfile, bHaveT,
                              opt2parg_bSet("-b", NPA, pa), tb,
                              opt2parg_bSet("-e", NPA, pa), te,
                              nreplica, &nset_d, &n_d, &dt_d, &remd.time);
    printf("Read %d sets of %d points in %s, dt = %g\n\n", nset_d, n_d, dfile, dt_d);

    if ((nset_t != nset_d) || (n_t != n_d) || (dt_t != dt_d))
    {
        gmx_fatal(FARGS, "Files %s and %s are inconsistent. Check log file",
                  tfile, dfile);
    }

    if (dfile2)
    {
        remd.data2 = read_xvg_time(dfile2, bHaveT,
                                   opt2parg_bSet("-b", NPA, pa), tb,
                                   opt2parg_bSet("-e", NPA, pa), te,
                                   nreplica, &nset_d2, &n_d2, &dt_d2, &remd.time);
        printf("Read %d sets of %d points in %s, dt = %g\n\n",
               nset_d2, n_d2, dfile2, dt_d2);
        if ((nset_d2 != nset_d) || (n_d != n_d2) || (dt_d != dt_d2))
        {
            gmx_fatal(FARGS, "Files %s and %s are inconsistent. Check log file",
                      dfile, dfile2);
        }
    }
    else
    {
        remd.data2 = NULL;
    }

    remd.nreplica  = nset_d;
    remd.nframe    = n_d;
    remd.dt        = 1;
    preprocess_remd(fp, &remd, cutoff, tref, ucut, bBack, Euf, Efu, Ei, t0, t1,
                    bSum, bDiscrete, nmult);

    optimize_remd_parameters(fp, &remd, maxiter, tol);

    dump_remd_parameters(fp, &remd, opt2fn("-o", NFILE, fnm),
                         opt2fn_null("-o2", NFILE, fnm),
                         opt2fn_null("-o3", NFILE, fnm),
                         opt2fn_null("-ee", NFILE, fnm),
                         opt2fn("-m", NFILE, fnm), skip, tref, oenv);

    if (bSplit)
    {
        printf("Splitting set of replicas in two halves\n");
        for (i = 0; (i < remd.nreplica); i++)
        {
            remd.bMask[i] = FALSE;
        }
        remd.nmask = 0;
        for (i = 0; (i < remd.nreplica); i += 2)
        {
            remd.bMask[i] = TRUE;
            remd.nmask++;
        }
        sum_ft(&remd);
        optimize_remd_parameters(fp, &remd, maxiter, tol);
        dump_remd_parameters(fp, &remd, "test1.xvg", NULL, NULL, NULL, NULL, skip, tref, oenv);

        for (i = 0; (i < remd.nreplica); i++)
        {
            remd.bMask[i] = !remd.bMask[i];
        }
        remd.nmask = remd.nreplica - remd.nmask;

        sum_ft(&remd);
        optimize_remd_parameters(fp, &remd, maxiter, tol);
        dump_remd_parameters(fp, &remd, "test2.xvg", NULL, NULL, NULL, NULL, skip, tref, oenv);

        for (i = 0; (i < remd.nreplica); i++)
        {
            remd.bMask[i] = FALSE;
        }
        remd.nmask = 0;
        for (i = 0; (i < remd.nreplica/2); i++)
        {
            remd.bMask[i] = TRUE;
            remd.nmask++;
        }
        sum_ft(&remd);
        optimize_remd_parameters(fp, &remd, maxiter, tol);
        dump_remd_parameters(fp, &remd, "test1.xvg", NULL, NULL, NULL, NULL, skip, tref, oenv);

        for (i = 0; (i < remd.nreplica); i++)
        {
            remd.bMask[i] = FALSE;
        }
        remd.nmask = 0;
        for (i = remd.nreplica/2; (i < remd.nreplica); i++)
        {
            remd.bMask[i] = TRUE;
            remd.nmask++;
        }
        sum_ft(&remd);
        optimize_remd_parameters(fp, &remd, maxiter, tol);
        dump_remd_parameters(fp, &remd, "test1.xvg", NULL, NULL, NULL, NULL, skip, tref, oenv);
    }
    ffclose(fp);

    view_all(oenv, NFILE, fnm);

    thanx(stderr);
#else
    fprintf(stderr, "This program should be compiled with the GNU scientific library. Please install the library and reinstall GROMACS.\n");
#endif /*HAVE_LIBGSL*/

    return 0;
}
