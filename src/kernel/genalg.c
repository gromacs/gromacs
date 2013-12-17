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
/**H*O*C**************************************************************
**                                                                  **
** No.!Version! Date ! Request !    Modification           ! Author **
** ---+-------+------+---------+---------------------------+------- **
**  1 + 3.1  +5/18/95+   -     + strategy DE/rand-to-best/1+  Storn **
**    +      +       +         + included                  +        **
**  1 + 3.2  +6/06/95+C.Fleiner+ change loops into memcpy  +  Storn **
**  2 + 3.2  +6/06/95+   -     + update comments           +  Storn **
**  1 + 3.3  +6/15/95+ K.Price + strategy DE/best/2 incl.  +  Storn **
**  2 + 3.3  +6/16/95+   -     + comments and beautifying  +  Storn **
**  3 + 3.3  +7/13/95+   -     + upper and lower bound for +  Storn **
**    +      +       +         + initialization            +        **
**  1 + 3.4  +2/12/96+   -     + increased printout prec.  +  Storn **
**  1 + 3.5  +5/28/96+   -     + strategies revisited      +  Storn **
**  2 + 3.5  +5/28/96+   -     + strategy DE/rand/2 incl.  +  Storn **
**  1 + 3.6  +8/06/96+ K.Price + Binomial Crossover added  +  Storn **
**  2 + 3.6  +9/30/96+ K.Price + cost variance output      +  Storn **
**  3 + 3.6  +9/30/96+   -     + alternative to ASSIGND    +  Storn **
**  4 + 3.6  +10/1/96+   -    + variable checking inserted +  Storn **
**  5 + 3.6  +10/1/96+   -     + strategy indic. improved  +  Storn **
**                                                                  **
***H*O*C*E***********************************************************/

/* Adopted for use in GROMACS by David van der Spoel, Oct 2001 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#include "typedefs.h"
#include "smalloc.h"
#include "futil.h"
#include "genalg.h"
#include "gmx_fatal.h"
#include "random.h"
#include "txtdump.h"
#include "vec.h"
#include "main.h"
#include "gmxfio.h"

static const char  *strat[] = {
    "DE/best/1/exp",          "DE/rand/1/exp",
    "DE/rand-to-best/1/exp",  "DE/best/2/exp",
    "DE/rand/2/exp",          "DE/best/1/bin",
    "DE/rand/1/bin",          "DE/rand-to-best/1/bin",
    "DE/best/2/bin",          "DE/rand/2/bin"
};

/*------------------------Macros----------------------------------------*/
static void assignd(int D, real a[], real b[])
{
    int i;

    for (i = 0; (i < D); i++)
    {
        a[i] = b[i];
    }
}

static real **make2d(int n, int m)
{
    int    i;
    real **r;

    snew(r, n);
    for (i = 0; (i < n); i++)
    {
        snew(r[i], m);
    }
    return r;
}

static void copy2range(int D, real c[], t_range r[])
{
    int i;

    for (i = 0; (i < D); i++)
    {
        /* Range check */
        while (c[i] < r[i].rmin)
        {
            c[i] += r[i].dr;
        }
        while (c[i] > r[i].rmax)
        {
            c[i] -= r[i].dr;
        }
        /*    if (c[i] < r[i].rmin)
           c[i] = r[i].rmin;
           if (c[i] > r[i].rmax)
           c[i] = r[i].rmax;
         */
        r[i].rval = c[i];
    }
}

t_genalg *init_ga(FILE *fplog, const char *infile, int D, t_range range[])
{
    FILE     *fpin_ptr;
    t_genalg *ga;
    double    ff = 0, cr = 0;
    int       i, j;

    /*------Initializations----------------------------*/
    snew(ga, 1);
    /*-----Read input data------------------------------------------------*/
    fpin_ptr   = gmx_fio_fopen(infile, "r");
    if (fscanf(fpin_ptr, "%d", &ga->NP) != 1         ||      /*---choice of strategy---*/
        fscanf(fpin_ptr, "%d", &ga->strategy) != 1   ||      /*---choice of strategy---*/
        fscanf(fpin_ptr, "%lf", &ff) != 1            ||      /*---weight factor------------*/
        fscanf(fpin_ptr, "%lf", &cr) != 1            ||      /*---crossing over factor-----*/
        fscanf(fpin_ptr, "%d", &ga->seed) != 1       ||      /*---random seed----------*/
        gmx_fio_fclose(fpin_ptr) != 1)
    {
        gmx_fatal(FARGS, "Error reading from file %s", infile);
    }

    ga->FF   = ff;
    ga->CR   = cr;
    ga->D    = D;
    ga->ipop = 0;
    ga->gen  = 0;

    /* Allocate memory */
    ga->pold = make2d(ga->NP, ga->D);
    ga->pnew = make2d(ga->NP, ga->D);
    snew(ga->tmp, ga->D);
    snew(ga->best, ga->D);
    snew(ga->bestit, ga->D);
    snew(ga->cost, ga->NP);
    snew(ga->msf, ga->NP);
    snew(ga->pres, ga->NP);
    snew(ga->scale, ga->NP);
    snew(ga->energy, ga->NP);

    /*-----Checking input variables for proper range--------------*/
    if ((ga->CR < 0) || (ga->CR > 1.0))
    {
        gmx_fatal(FARGS, "CR=%f, should be ex [0,1]", ga->CR);
    }
    if (ga->seed <= 0)
    {
        gmx_fatal(FARGS, "seed=%d, should be > 0", ga->seed);
    }
    if ((ga->strategy < 0) || (ga->strategy > 10))
    {
        gmx_fatal(FARGS, "strategy=%d, should be ex {1-10}", ga->strategy);
    }

    /* spread initial population members */
    for (i = 0; (i < ga->NP); i++)
    {
        for (j = 0; (j < ga->D); j++)
        {
            ga->pold[i][j] = value_rand(&(range[j]), &ga->seed);
        }
    }

    fprintf(fplog, "-----------------------------------------------\n");
    fprintf(fplog, "Genetic algorithm parameters\n");
    fprintf(fplog, "-----------------------------------------------\n");
    fprintf(fplog, "Number of variables:   %d\n", ga->D);
    fprintf(fplog, "Population size:       %d\n", ga->NP);
    fprintf(fplog, "Strategy:              %s\n", strat[ga->strategy]);
    fprintf(fplog, "Weight factor:         %g\n", ga->FF);
    fprintf(fplog, "Crossing over factor:  %g\n", ga->CR);
    fprintf(fplog, "Random seed:           %d\n", ga->seed);
    fprintf(fplog, "-----------------------------------------------\n");

    return ga;
}

void update_ga(FILE *fpout_ptr, t_range range[], t_genalg *ga)
{
    static int  i_init = 0;         /* Initialisation related stuff       */
    int         i, j, L, n;         /* counting variables                 */
    int         r1, r2, r3, r4, r5; /* placeholders for random indexes    */

    if (i_init < ga->NP)
    {
        /* Copy data for first force evaluation to range array  */
        copy2range(ga->D, ga->pold[i_init], range);

        i_init++;
        return;
    }
    else
    {
        /* Now starts real genetic stuff, a new trial set is made */
        if (ga->ipop == ga->NP)
        {
            ga->gen++;
            i = ga->ipop = 0;
        }
        else
        {
            i = ga->ipop;
        }

        do /* Pick a random population member */
        {  /* Endless loop for ga->NP < 2 !!!     */
            r1 = (int)(rando(&ga->seed)*ga->NP);
        }
        while (r1 == i);

        do /* Pick a random population member */
        {  /* Endless loop for ga->NP < 3 !!!     */
            r2 = (int)(rando(&ga->seed)*ga->NP);
        }
        while ((r2 == i) || (r2 == r1));

        do
        {
            /* Pick a random population member */
            /* Endless loop for ga->NP < 4 !!!     */
            r3 = (int)(rando(&ga->seed)*ga->NP);
        }
        while ((r3 == i) || (r3 == r1) || (r3 == r2));

        do
        {
            /* Pick a random population member */
            /* Endless loop for ga->NP < 5 !!!     */
            r4 = (int)(rando(&ga->seed)*ga->NP);
        }
        while ((r4 == i) || (r4 == r1) || (r4 == r2) || (r4 == r3));

        do
        {
            /* Pick a random population member */
            /* Endless loop for ga->NP < 6 !!!     */
            r5 = (int)(rando(&ga->seed)*ga->NP);
        }
        while ((r5 == i) || (r5 == r1) || (r5 == r2) || (r5 == r3) || (r5 == r4));


        /* Choice of strategy
         * We have tried to come up with a sensible naming-convention: DE/x/y/z
         * DE :  stands for Differential Evolution
         * x  :  a string which denotes the vector to be perturbed
         * y  :  number of difference vectors taken for perturbation of x
         * z  :  crossover method (exp = exponential, bin = binomial)
         *
         * There are some simple rules which are worth following:
         * 1)  ga->FF is usually between 0.5 and 1 (in rare cases > 1)
         * 2)  ga->CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to
         *     be tried first
         * 3)  To start off ga->NP = 10*ga->D is a reasonable choice. Increase ga->NP if
         *     misconvergence happens.
         * 4)  If you increase ga->NP, ga->FF usually has to be decreased
         * 5)  When the DE/ga->best... schemes fail DE/rand... usually works and
         *     vice versa
         * EXPONENTIAL ga->CROSSOVER
         *-------DE/ga->best/1/exp-------
         *-------Our oldest strategy but still not bad. However, we have found several
         *-------optimization problems where misconvergence occurs.
         */
        assignd(ga->D, ga->tmp, ga->pold[i]);
        n = (int)(rando(&ga->seed)*ga->D);
        L = 0;

        switch (ga->strategy)
        {
            case 1:
                /* strategy DE0 (not in our paper) */
                do
                {
                    ga->tmp[n] = ga->bestit[n] + ga->FF*(ga->pold[r2][n]-ga->pold[r3][n]);
                    n          = (n+1)%ga->D;
                    L++;
                }
                while ((rando(&ga->seed) < ga->CR) && (L < ga->D));
                break;

            /* DE/rand/1/exp
             * This is one of my favourite strategies. It works especially
             * well when the ga->bestit[]"-schemes experience misconvergence.
             * Try e.g. ga->FF=0.7 and ga->CR=0.5 * as a first guess.
             */
            case 2:
                /* strategy DE1 in the techreport */
                do
                {
                    ga->tmp[n] = ga->pold[r1][n] + ga->FF*(ga->pold[r2][n]-ga->pold[r3][n]);
                    n          = (n+1)%ga->D;
                    L++;
                }
                while ((rando(&ga->seed) < ga->CR) && (L < ga->D));
                break;

            /* DE/rand-to-ga->best/1/exp
             * This strategy seems to be one of the ga->best strategies.
             * Try ga->FF=0.85 and ga->CR=1.
             * If you get misconvergence try to increase ga->NP.
             * If this doesn't help you should play around with all three
             * control variables.
             */
            case 3:
                /* similar to DE2 but generally better */
                do
                {
                    ga->tmp[n] = ga->tmp[n] + ga->FF*(ga->bestit[n] - ga->tmp[n]) +
                        ga->FF*(ga->pold[r1][n]-ga->pold[r2][n]);
                    n = (n+1)%ga->D;
                    L++;
                }
                while ((rando(&ga->seed) < ga->CR) && (L < ga->D));
                break;

            /* DE/ga->best/2/exp is another powerful strategy worth trying */
            case 4:
                do
                {
                    ga->tmp[n] = ga->bestit[n] +
                        (ga->pold[r1][n]+ga->pold[r2][n]-ga->pold[r3][n]-ga->pold[r4][n])*ga->FF;
                    n = (n+1)%ga->D;
                    L++;
                }
                while ((rando(&ga->seed) < ga->CR) && (L < ga->D));
                break;

            /*----DE/rand/2/exp seems to be a robust optimizer for many functions-----*/
            case 5:
                do
                {
                    ga->tmp[n] = ga->pold[r5][n] +
                        (ga->pold[r1][n]+ga->pold[r2][n]-ga->pold[r3][n]-ga->pold[r4][n])*ga->FF;
                    n = (n+1)%ga->D;
                    L++;
                }
                while ((rando(&ga->seed) < ga->CR) && (L < ga->D));
                break;

            /*===Essentially same strategies but BINOMIAL ga->CROSSOVER===*/

            /*-------DE/ga->best/1/bin------*/
            case 6:
                for (L = 0; L < ga->D; L++)
                {
                    /* perform D binomial trials */
                    if ((rando(&ga->seed) < ga->CR) || (L == (ga->D-1)))
                    {
                        /* change at least one parameter */
                        ga->tmp[n] = ga->bestit[n] + ga->FF*(ga->pold[r2][n]-ga->pold[r3][n]);
                    }
                    n = (n+1)%ga->D;
                }
                break;

            /*-------DE/rand/1/bin------*/
            case 7:
                for (L = 0; L < ga->D; L++)
                {
                    /* perform D binomial trials */
                    if ((rando(&ga->seed) < ga->CR) || (L == (ga->D-1)))
                    {
                        /* change at least one parameter */
                        ga->tmp[n] = ga->pold[r1][n] + ga->FF*(ga->pold[r2][n]-ga->pold[r3][n]);
                    }
                    n = (n+1)%ga->D;
                }
                break;

            /*-------DE/rand-to-ga->best/1/bin------*/
            case 8:
                for (L = 0; L < ga->D; L++)
                {
                    /* perform ga->D binomial trials */
                    if ((rando(&ga->seed) < ga->CR) || (L == (ga->D-1)))
                    {
                        /* change at least one parameter */
                        ga->tmp[n] = ga->tmp[n] + ga->FF*(ga->bestit[n] - ga->tmp[n]) +
                            ga->FF*(ga->pold[r1][n]-ga->pold[r2][n]);
                    }
                    n = (n+1)%ga->D;
                }
                break;

            /*-------DE/ga->best/2/bin------*/
            case 9:
                for (L = 0; L < ga->D; L++)
                {
                    /* perform ga->D binomial trials */
                    if ((rando(&ga->seed) < ga->CR) || (L == (ga->D-1)))
                    {
                        /* change at least one parameter */
                        ga->tmp[n] = ga->bestit[n] +
                            (ga->pold[r1][n]+ga->pold[r2][n]-ga->pold[r3][n]-ga->pold[r4][n])*ga->FF;
                    }
                    n = (n+1)%ga->D;
                }
                break;

            /*-------DE/rand/2/bin-------*/
            default:
                for (L = 0; L < ga->D; L++)
                {
                    /* perform ga->D binomial trials */
                    if ((rando(&ga->seed) < ga->CR) || (L == (ga->D-1)))
                    {
                        /* change at least one parameter */
                        ga->tmp[n] = ga->pold[r5][n] +
                            (ga->pold[r1][n]+ga->pold[r2][n]-ga->pold[r3][n]-ga->pold[r4][n])*ga->FF;
                    }
                    n = (n+1)%ga->D;
                }
                break;
        }

        /*===Trial mutation now in ga->tmp[]. Test how good this choice really was.==*/
        copy2range(ga->D, ga->tmp, range);
    }
}

gmx_bool print_ga(FILE *fp, t_genalg *ga, real msf, tensor pres, rvec scale,
                  real energy, t_range range[], real tol)
{
    static int      nfeval = 0; /* number of function evaluations     */
    static gmx_bool bImproved;
    real            trial_cost;
    real            cvar;  /* computes the cost variance         */
    real            cmean; /* mean cost                          */
    int             i, j;
    real          **pswap;

    trial_cost = cost(pres, msf, energy);
    if (nfeval < ga->NP)
    {
        ga->cost[nfeval]   = trial_cost;
        ga->msf[nfeval]    = msf;
        ga->energy[nfeval] = energy;
        copy_mat(pres, ga->pres[nfeval]);
        copy_rvec(scale, ga->scale[nfeval]);
        if (debug)
        {
            pr_rvec(debug, 0, "scale", scale, DIM, TRUE);
            pr_rvec(debug, 0, "pold ", ga->pold[nfeval]+4, DIM, TRUE);
        }
        nfeval++;
        return FALSE;
    }
    /* When we get here we have done an initial evaluation for all
     * animals in the population
     */
    if (ga->ipop == 0)
    {
        bImproved = FALSE;
    }

    /* First iteration after first round of trials */
    if (nfeval == ga->NP)
    {
        /* Evaluate who is ga->best */
        ga->imin = 0;
        for (j = 1; (j < ga->NP); j++)
        {
            if (ga->cost[j] < ga->cost[ga->imin])
            {
                ga->imin = j;
            }
        }
        assignd(ga->D, ga->best, ga->pold[ga->imin]);
        /* save best member ever          */
        assignd(ga->D, ga->bestit, ga->pold[ga->imin]);
        /* save best member of generation */
    }

    if (trial_cost < ga->cost[ga->ipop])
    {
        if (trial_cost < ga->cost[ga->imin])
        {
            /* Was this a new minimum? */
            /* if so, reset cmin to new low...*/
            ga->imin = ga->ipop;
            assignd(ga->D, ga->best, ga->tmp);
            bImproved = TRUE;
        }
        /* improved objective function value ? */
        ga->cost[ga->ipop]   = trial_cost;

        ga->msf[ga->ipop]    = msf;
        ga->energy[ga->ipop] = energy;
        copy_mat(pres, ga->pres[ga->ipop]);
        copy_rvec(scale, ga->scale[ga->ipop]);

        assignd(ga->D, ga->pnew[ga->ipop], ga->tmp);

    }
    else
    {
        /* replace target with old value */
        assignd(ga->D, ga->pnew[ga->ipop], ga->pold[ga->ipop]);
    }
    /* #define SCALE_DEBUG*/
#ifdef SCALE_DEBUG
    if (ga->D > 5)
    {
        rvec dscale;
        rvec_sub(ga->scale[ga->imin], &(ga->best[ga->D-3]), dscale);
        if (norm(dscale) > 0)
        {
            pr_rvec(fp, 0, "scale", scale, DIM, TRUE);
            pr_rvec(fp, 0, "best ", &(ga->best[ga->D-3]), DIM, TRUE);
            fprintf(fp, "imin = %d, ipop = %d, nfeval = %d\n", ga->imin,
                    ga->ipop, nfeval);
            gmx_fatal(FARGS, "Scale inconsistency");
        }
    }
#endif

    /* Increase population member count */
    ga->ipop++;

    /* End mutation loop through population? */
    if (ga->ipop == ga->NP)
    {
        /* Save ga->best population member of current iteration */
        assignd(ga->D, ga->bestit, ga->best);

        /* swap population arrays. New generation becomes old one */
        pswap     = ga->pold;
        ga->pold  = ga->pnew;
        ga->pnew  = pswap;

        /*----Compute the energy variance (just for monitoring purposes)-----------*/
        /* compute the mean value first */
        cmean = 0.0;
        for (j = 0; (j < ga->NP); j++)
        {
            cmean += ga->cost[j];
        }
        cmean = cmean/ga->NP;

        /* now the variance              */
        cvar = 0.0;
        for (j = 0; (j < ga->NP); j++)
        {
            cvar += sqr(ga->cost[j] - cmean);
        }
        cvar = cvar/(ga->NP-1);

        /*----Output part----------------------------------------------------------*/
        if (1 || bImproved || (nfeval == ga->NP))
        {
            fprintf(fp, "\nGen: %5d\n  Cost: %12.3e  <Cost>: %12.3e\n"
                    "  Ener: %8.4f RMSF: %8.3f Pres: %8.1f %8.1f  %8.1f (%8.1f)\n"
                    "  Box-Scale: %15.10f  %15.10f  %15.10f\n",
                    ga->gen, ga->cost[ga->imin], cmean, ga->energy[ga->imin],
                    sqrt(ga->msf[ga->imin]), ga->pres[ga->imin][XX][XX],
                    ga->pres[ga->imin][YY][YY], ga->pres[ga->imin][ZZ][ZZ],
                    trace(ga->pres[ga->imin])/3.0,
                    ga->scale[ga->imin][XX], ga->scale[ga->imin][YY],
                    ga->scale[ga->imin][ZZ]);

            for (j = 0; (j < ga->D); j++)
            {
                fprintf(fp, "\tbest[%d]=%-15.10f\n", j, ga->best[j]);
            }
            if (ga->cost[ga->imin] < tol)
            {
                for (i = 0; (i < ga->NP); i++)
                {
                    fprintf(fp, "Animal: %3d Cost:%8.3f  Ener: %8.3f RMSF: %8.3f%s\n",
                            i, ga->cost[i], ga->energy[i], sqrt(ga->msf[i]),
                            (i == ga->imin) ? " ***" : "");
                    for (j = 0; (j < ga->D); j++)
                    {
                        fprintf(fp, "\tParam[%d][%d]=%15.10g\n", i, j, ga->pold[i][j]);
                    }
                }
                return TRUE;
            }
            fflush(fp);
        }
    }
    nfeval++;

    return FALSE;
}
