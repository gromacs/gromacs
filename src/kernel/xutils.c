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

#include <assert.h>

#include "typedefs.h"
#include "smalloc.h"
#include "strdb.h"
#include "string2.h"
#include "xmdrun.h"
#include "vec.h"
#include "genalg.h"
#include "random.h"

real mol_dipole(int k0, int k1, rvec x[], real q[])
{
    int  k, m;
    rvec mu;

    clear_rvec(mu);
    for (k = k0; (k < k1); k++)
    {
        for (m = 0; (m < DIM); m++)
        {
            mu[m] += q[k]*x[k][m];
        }
    }
    return norm(mu); /* Dipole moment of this molecule in e nm */
}

real calc_mu_aver(t_commrec *cr, rvec x[], real q[], rvec mu,
                  t_block *mols, t_mdatoms *md, int gnx, atom_id grpindex[])
{
    int     i, start, end;
    real    mu_ave;

    start = md->start;
    end   = md->homenr + start;

    /*
       clear_rvec(mu);
       for(i=start; (i<end); i++)
       for(m=0; (m<DIM); m++)
        mu[m] += q[i]*x[i][m];
       if (PAR(cr)) {
       gmx_sum(DIM,mu,cr);
       }
     */
    /* I guess we have to parallelise this one! */

    if (gnx > 0)
    {
        mu_ave = 0.0;
        for (i = 0; (i < gnx); i++)
        {
            int gi = grpindex[i];
            mu_ave += mol_dipole(mols->index[gi], mols->index[gi+1], x, q);
        }

        return(mu_ave/gnx);
    }
    else
    {
        return 0;
    }
}

/* Lots of global variables! Yummy... */
static t_ffscan ff;

void set_ffvars(t_ffscan *fff)
{
    ff = *fff;
}

real cost(tensor P, real MSF, real E)
{
    return (ff.fac_msf*MSF+ff.fac_epot*sqr(E-ff.epot)+ff.fac_pres*
            (sqr(P[XX][XX]-ff.pres)+sqr(P[YY][YY]-ff.pres)+sqr(P[ZZ][ZZ]-ff.pres)));

}

static const char     *esenm[eseNR] = { "SIG", "EPS", "BHAMA", "BHAMB", "BHAMC", "CELlX", "CELLY", "CELLZ" };
static int             nparm = 0, *param_val = NULL;
static t_range        *range = NULL;
static t_genalg       *ga    = NULL;
static rvec            scale = { 1, 1, 1 };

static void init_range(t_range *r, int np, int atype, int ptype,
                       real rmin, real rmax)
{
    if (rmin > rmax)
    {
        gmx_fatal(FARGS, "rmin (%f) > rmax (%f)", rmin, rmax);
    }
    if (np <= 0)
    {
        gmx_fatal(FARGS, "np (%d) should be > 0", np);
    }
    if ((rmax > rmin) && (np <= 1))
    {
        gmx_fatal(FARGS, "If rmax > rmin, np should be > 1");
    }
    if ((ptype < 0) || (ptype >= eseNR))
    {
        gmx_fatal(FARGS, "ptype (%d) should be < %d", ptype, eseNR);
    }
    r->np    = np;
    r->atype = atype;
    r->ptype = ptype;
    r->rmin  = rmin;
    r->rmax  = rmax;
    r->rval  = rmin;
    r->dr    = r->rmax - r->rmin;
}

static t_range *read_range(const char *db, int *nrange)
{
    int       nlines, nr, np, i;
    char    **lines;
    t_range  *range;
    int       atype, ptype;
    double    rmin, rmax;

    nlines = get_file(db, &lines);
    snew(range, nlines);

    nr = 0;
    for (i = 0; (i < nlines); i++)
    {
        strip_comment(lines[i]);
        if (sscanf(lines[i], "%d%d%d%lf%lf", &np, &atype, &ptype, &rmin, &rmax) == 5)
        {
            if (ff.bLogEps && (ptype == eseEPSILON) && (rmin <= 0))
            {
                gmx_fatal(FARGS, "When using logarithmic epsilon increments the minimum"
                          "value must be > 0");
            }
            init_range(&range[nr], np, atype, ptype, rmin, rmax);
            nr++;
        }
    }
    fprintf(stderr, "found %d variables to iterate over\n", nr);

    *nrange = nr;

    for (nr = 0; (nr < nlines); nr++)
    {
        sfree(lines[nr]);
    }
    sfree(lines);

    return range;
}

static real value_range(t_range *r, int n)
{
    real logrmin, logrmax;

    if ((n < 0) || (n > r->np))
    {
        gmx_fatal(FARGS, "Value (%d) out of range for value_range (max %d)", n, r->np);
    }

    if (r->np == 1)
    {
        r->rval = r->rmin;
    }
    else
    {
        if ((r->ptype == eseEPSILON) && ff.bLogEps)
        {
            logrmin = log(r->rmin);
            logrmax = log(r->rmax);
            r->rval = exp(logrmin + (n*(logrmax-logrmin))/(r->np-1));
        }
        else
        {
            r->rval = r->rmin+(n*(r->dr))/(r->np-1);
        }
    }
    return r->rval;
}

real value_rand(t_range *r, int *seed)
{
    real logrmin, logrmax;
    real mr;

    if (r->np == 1)
    {
        r->rval = r->rmin;
    }
    else
    {
        mr = rando(seed);
        if ((r->ptype == eseEPSILON) && ff.bLogEps)
        {
            logrmin = log(r->rmin);
            logrmax = log(r->rmax);
            r->rval = exp(logrmin + mr*(logrmax-logrmin));
        }
        else
        {
            r->rval = r->rmin + mr*(r->rmax-r->rmin);
        }
    }
    if (debug)
    {
        fprintf(debug, "type: %s, value: %g\n", esenm[r->ptype], r->rval);
    }
    return r->rval;
}

static void update_ff(t_forcerec *fr, int nparm, t_range range[], int param_val[])
{
    static double *sigma = NULL, *eps = NULL, *c6 = NULL, *cn = NULL, *bhama = NULL, *bhamb = NULL, *bhamc = NULL;
    real           val, *nbfp;
    int            i, j, atnr;

    atnr = fr->ntype;
    nbfp = fr->nbfp;

    if (fr->bBHAM)
    {
        if (bhama == NULL)
        {
            assert(bhamb == NULL && bhamc == NULL);
            snew(bhama, atnr);
            snew(bhamb, atnr);
            snew(bhamc, atnr);
        }
    }
    else
    {
        if (sigma == NULL)
        {
            assert(eps == NULL && c6 == NULL && cn == NULL);
            snew(sigma, atnr);
            snew(eps, atnr);
            snew(c6, atnr);
            snew(cn, atnr);
        }
    }
    /* Get current values for everything */
    for (i = 0; (i < nparm); i++)
    {
        if (ga)
        {
            val = range[i].rval;
        }
        else
        {
            val = value_range(&range[i], param_val[i]);
        }
        if (debug)
        {
            fprintf(debug, "val = %g\n", val);
        }
        switch (range[i].ptype)
        {
            case eseSIGMA:
                assert(!fr->bBHAM);
                sigma[range[i].atype] = val;
                break;
            case eseEPSILON:
                assert(!fr->bBHAM);
                eps[range[i].atype] = val;
                break;
            case eseBHAMA:
                assert(fr->bBHAM);
                bhama[range[i].atype] = val;
                break;
            case eseBHAMB:
                assert(fr->bBHAM);
                bhamb[range[i].atype] = val;
                break;
            case eseBHAMC:
                assert(fr->bBHAM);
                bhamc[range[i].atype] = val;
                break;
            case eseCELLX:
                scale[XX] = val;
                break;
            case eseCELLY:
                scale[YY] = val;
                break;
            case eseCELLZ:
                scale[ZZ] = val;
                break;
            default:
                gmx_fatal(FARGS, "Unknown ptype");
        }
    }
    if (fr->bBHAM)
    {
        for (i = 0; (i < atnr); i++)
        {
            for (j = 0; (j <= i); j++)
            {
                BHAMA(nbfp, atnr, i, j) = BHAMA(nbfp, atnr, j, i) = sqrt(bhama[i]*bhama[j]);
                BHAMB(nbfp, atnr, i, j) = BHAMB(nbfp, atnr, j, i) = sqrt(bhamb[i]*bhamb[j]);
                BHAMC(nbfp, atnr, i, j) = BHAMC(nbfp, atnr, j, i) = sqrt(bhamc[i]*bhamc[j]);
            }
        }
    }
    else
    {
        /* Now build a new matrix */
        for (i = 0; (i < atnr); i++)
        {
            c6[i] = 4*eps[i]*pow(sigma[i], 6.0);
            cn[i] = 4*eps[i]*pow(sigma[i], ff.npow);
        }
        for (i = 0; (i < atnr); i++)
        {
            for (j = 0; (j <= i); j++)
            {
                C6(nbfp, atnr, i, j)  = C6(nbfp, atnr, j, i)  = sqrt(c6[i]*c6[j]);
                C12(nbfp, atnr, i, j) = C12(nbfp, atnr, j, i) = sqrt(cn[i]*cn[j]);
            }
        }
    }

    if (debug)
    {
        if (!fr->bBHAM)
        {
            for (i = 0; (i < atnr); i++)
            {
                fprintf(debug, "atnr = %2d  sigma = %8.4f  eps = %8.4f\n", i, sigma[i], eps[i]);
            }
        }
        for (i = 0; (i < atnr); i++)
        {
            for (j = 0; (j < atnr); j++)
            {
                if (fr->bBHAM)
                {
                    fprintf(debug, "i: %2d  j: %2d  A:  %10.5e  B:  %10.5e  C:  %10.5e\n", i, j,
                            BHAMA(nbfp, atnr, i, j), BHAMB(nbfp, atnr, i, j), BHAMC(nbfp, atnr, i, j));
                }
                else
                {
                    fprintf(debug, "i: %2d  j: %2d  c6:  %10.5e  cn:  %10.5e\n", i, j,
                            C6(nbfp, atnr, i, j), C12(nbfp, atnr, i, j));
                }
            }
        }
    }
}

static void scale_box(int natoms, rvec x[], matrix box)
{
    int i, m;

    if ((scale[XX] != 1.0) ||   (scale[YY] != 1.0) ||   (scale[ZZ] != 1.0))
    {
        if (debug)
        {
            fprintf(debug, "scale = %8.4f  %8.4f  %8.4f\n",
                    scale[XX], scale[YY], scale[ZZ]);
        }
        for (m = 0; (m < DIM); m++)
        {
            box[m][m] *= scale[m];
        }
        for (i = 0; (i < natoms); i++)
        {
            for (m = 0; (m < DIM); m++)
            {
                x[i][m] *= scale[m];
            }
        }
    }
}

gmx_bool update_forcefield(FILE *fplog,
                           int nfile, const t_filenm fnm[], t_forcerec *fr,
                           int natoms, rvec x[], matrix box)
{
    static int ntry, ntried;
    int        i, j;
    gmx_bool   bDone;

    /* First time around we have to read the parameters */
    if (nparm == 0)
    {
        range = read_range(ftp2fn(efDAT, nfile, fnm), &nparm);
        if (nparm == 0)
        {
            gmx_fatal(FARGS, "No correct parameter info in %s", ftp2fn(efDAT, nfile, fnm));
        }
        snew(param_val, nparm);

        if (opt2bSet("-ga", nfile, fnm))
        {
            /* Genetic algorithm time */
            ga = init_ga(fplog, opt2fn("-ga", nfile, fnm), nparm, range);
        }
        else
        {
            /* Determine the grid size */
            ntry = 1;
            for (i = 0; (i < nparm); i++)
            {
                ntry *= range[i].np;
            }
            ntried = 0;

            fprintf(fplog, "Going to try %d different combinations of %d parameters\n",
                    ntry, nparm);
        }
    }
    if (ga)
    {
        update_ga(fplog, range, ga);
    }
    else
    {
        /* Increment the counter
         * Non-trivial, since this is nparm nested loops in principle
         */
        for (i = 0; (i < nparm); i++)
        {
            if (param_val[i] < (range[i].np-1))
            {
                param_val[i]++;
                for (j = 0; (j < i); j++)
                {
                    param_val[j] = 0;
                }
                ntried++;
                break;
            }
        }
        if (i == nparm)
        {
            fprintf(fplog, "Finished with %d out of %d iterations\n", ntried+1, ntry);
            return TRUE;
        }
    }

    /* Now do the real updating */
    update_ff(fr, nparm, range, param_val);

    /* Update box and coordinates if necessary */
    scale_box(natoms, x, box);

    return FALSE;
}

static void print_range(FILE *fp, tensor P, real MSF, real energy)
{
    int  i;

    fprintf(fp, "%8.3f  %8.3f  %8.3f  %8.3f",
            cost(P, MSF, energy), trace(P)/3, MSF, energy);
    for (i = 0; (i < nparm); i++)
    {
        fprintf(fp, " %s %10g", esenm[range[i].ptype], range[i].rval);
    }
    fprintf(fp, " FF\n");
    fflush(fp);
}

static real msf(int n, rvec f1[], rvec f2[])
{
    int  i, j;
    rvec ff2;
    real msf1 = 0;

    for (i = 0; (i < n); )
    {
        clear_rvec(ff2);
        for (j = 0; ((j < ff.molsize) && (i < n)); j++, i++)
        {
            rvec_inc(ff2, f1[i]);
            if (f2)
            {
                rvec_inc(ff2, f2[i]);
            }
        }
        msf1 += iprod(ff2, ff2);
    }

    return msf1/n;
}

static void print_grid(FILE *fp, real ener[], int natoms, rvec f[], rvec fshake[],
                       rvec x[], t_block *mols, real mass[], tensor pres)
{
    static gmx_bool    bFirst = TRUE;
    static const char *desc[] = {
        "------------------------------------------------------------------------",
        "In the output from the forcefield scan we have the potential energy,",
        "then the root mean square force on the atoms, and finally the parameters",
        "in the order they appear in the input file.",
        "------------------------------------------------------------------------"
    };
    real               msf1;
    int                i;

    if (bFirst)
    {
        for (i = 0; (i < asize(desc)); i++)
        {
            fprintf(fp, "%s\n", desc[i]);
        }
        fflush(fp);
        bFirst = FALSE;
    }
    if ((ff.tol == 0) || (fabs(ener[F_EPOT]/ff.nmol-ff.epot) < ff.tol))
    {
        msf1 = msf(natoms, f, fshake);
        if ((ff.f_max == 0) || (msf1 < sqr(ff.f_max)))
        {
            print_range(fp, pres, msf1, ener[F_EPOT]/ff.nmol);
        }
    }
}

gmx_bool print_forcefield(FILE *fp, real ener[], int natoms, rvec f[], rvec fshake[],
                          rvec x[], t_block *mols, real mass[], tensor pres)
{
    real msf1;

    if (ga)
    {
        msf1 = msf(natoms, f, fshake);
        if (debug)
        {
            fprintf(fp, "Pressure: %12g, RMSF: %12g, Energy-Epot: %12g, cost: %12g\n",
                    ener[F_PRES], sqrt(msf1), ener[F_EPOT]/ff.nmol-ff.epot,
                    cost(pres, msf1, ener[F_EPOT]/ff.nmol));
        }
        if (print_ga(fp, ga, msf1, pres, scale, (ener[F_EPOT]/ff.nmol), range, ff.tol))
        {
            return TRUE;
        }
        fflush(fp);
    }
    else
    {
        print_grid(fp, ener, natoms, f, fshake, x, mols, mass, pres);
    }
    return FALSE;
}
