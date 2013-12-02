/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "typedefs.h"
#include "gmx_fatal.h"
#include "mdrun.h"
#include "md_support.h"
#include "md_logging.h"
#include "types/iteratedconstraints.h"

#ifdef GMX_DOUBLE
#define CONVERGEITER  0.000000001
#define CLOSE_ENOUGH  0.000001000
#else
#define CONVERGEITER  0.0001
#define CLOSE_ENOUGH  0.0050
#endif

/* we want to keep track of the close calls.  If there are too many, there might be some other issues.
   so we make sure that it's either less than some predetermined number, or if more than that number,
   only some small fraction of the total. */
#define MAX_NUMBER_CLOSE        50
#define FRACTION_CLOSE       0.001

/* maximum length of cyclic traps to check, emerging from limited numerical precision  */
#define CYCLEMAX            20

void gmx_iterate_init(gmx_iterate_t *iterate, gmx_bool bSetIterationActive)
{
    int i;

    iterate->iter_i           = 0;
    iterate->bIterationActive = bSetIterationActive;
    iterate->num_close        = 0;
    for (i = 0; i < MAXITERCONST+2; i++)
    {
        iterate->allrelerr[i] = 0;
    }
}

gmx_bool done_iterating(const t_commrec *cr, FILE *fplog, int nsteps, gmx_iterate_t *iterate, gmx_bool bFirstIterate, real fom, real *newf)
{
    /* monitor convergence, and use a secant search to propose new
       values.
                                                                  x_{i} - x_{i-1}
       The secant method computes x_{i+1} = x_{i} - f(x_{i}) * ---------------------
                                                                f(x_{i}) - f(x_{i-1})

       The function we are trying to zero is fom-x, where fom is the
       "figure of merit" which is the pressure (or the veta value) we
       would get by putting in an old value of the pressure or veta into
       the incrementor function for the step or half step.  I have
       verified that this gives the same answer as self consistent
       iteration, usually in many fewer steps, especially for small tau_p.

       We could possibly eliminate an iteration with proper use
       of the value from the previous step, but that would take a bit
       more bookkeeping, especially for veta, since tests indicate the
       function of veta on the last step is not sufficiently close to
       guarantee convergence this step. This is
       good enough for now.  On my tests, I could use tau_p down to
       0.02, which is smaller that would ever be necessary in
       practice. Generally, 3-5 iterations will be sufficient */

    real     relerr, err, xmin;
    int      i;
    gmx_bool incycle;

    if (bFirstIterate)
    {
        iterate->x     = fom;
        iterate->f     = fom-iterate->x;
        iterate->xprev = 0;
        iterate->fprev = 0;
        *newf          = fom;
    }
    else
    {
        iterate->f = fom-iterate->x; /* we want to zero this difference */
        if ((iterate->iter_i > 1) && (iterate->iter_i < MAXITERCONST))
        {
            if (iterate->f == iterate->fprev)
            {
                *newf = iterate->f;
            }
            else
            {
                *newf = iterate->x - (iterate->x-iterate->xprev)*(iterate->f)/(iterate->f-iterate->fprev);
            }
        }
        else
        {
            /* just use self-consistent iteration the first step to initialize, or
               if it's not converging (which happens occasionally -- need to investigate why) */
            *newf = fom;
        }
    }
    /* Consider a slight shortcut allowing us to exit one sooner -- we check the
       difference between the closest of x and xprev to the new
       value. To be 100% certain, we should check the difference between
       the last result, and the previous result, or

       relerr = (fabs((x-xprev)/fom));

       but this is pretty much never necessary under typical conditions.
       Checking numerically, it seems to lead to almost exactly the same
       trajectories, but there are small differences out a few decimal
       places in the pressure, and eventually in the v_eta, but it could
       save an interation.

       if (fabs(*newf-x) < fabs(*newf - xprev)) { xmin = x;} else { xmin = xprev;}
       relerr = (fabs((*newf-xmin) / *newf));
     */

    err    = fabs((iterate->f-iterate->fprev));
    relerr = fabs(err/fom);

    iterate->allrelerr[iterate->iter_i] = relerr;

    if (iterate->iter_i > 0)
    {
        if (debug)
        {
            fprintf(debug, "Iterating NPT constraints: %6i %20.12f%14.6g%20.12f\n",
                    iterate->iter_i, fom, relerr, *newf);
        }

        if ((relerr < CONVERGEITER) || (err < CONVERGEITER) || (fom == 0) || ((iterate->x == iterate->xprev) && iterate->iter_i > 1))
        {
            iterate->bIterationActive = FALSE;
            if (debug)
            {
                fprintf(debug, "Iterating NPT constraints: CONVERGED\n");
            }
            return TRUE;
        }
        if (iterate->iter_i > MAXITERCONST)
        {
            if (relerr < CLOSE_ENOUGH)
            {
                incycle = FALSE;
                for (i = 1; i < CYCLEMAX; i++)
                {
                    if ((iterate->allrelerr[iterate->iter_i-(1+i)] == iterate->allrelerr[iterate->iter_i-1]) &&
                        (iterate->allrelerr[iterate->iter_i-(1+i)] == iterate->allrelerr[iterate->iter_i-(1+2*i)]))
                    {
                        incycle = TRUE;
                        if (debug)
                        {
                            fprintf(debug, "Exiting from an NPT iterating cycle of length %d\n", i);
                        }
                        break;
                    }
                }

                if (incycle)
                {
                    /* step 1: trapped in a numerical attractor */
                    /* we are trapped in a numerical attractor, and can't converge any more, and are close to the final result.
                       Better to give up convergence here than have the simulation die.
                     */
                    iterate->num_close++;
                    iterate->bIterationActive = FALSE;
                    return TRUE;
                }
                else
                {
                    /* Step #2: test if we are reasonably close for other reasons, then monitor the number.  If not, die */

                    /* how many close calls have we had?  If less than a few, we're OK */
                    if (iterate->num_close < MAX_NUMBER_CLOSE)
                    {
                        md_print_warn(cr, fplog, "Slight numerical convergence deviation with NPT at step %d, relative error only %10.5g, likely not a problem, continuing\n", nsteps, relerr);
                        iterate->num_close++;
                        iterate->bIterationActive = FALSE;
                        return TRUE;
                        /* if more than a few, check the total fraction.  If too high, die. */
                    }
                    else if (iterate->num_close/(double)nsteps > FRACTION_CLOSE)
                    {
                        gmx_fatal(FARGS, "Could not converge NPT constraints, too many exceptions (%d%%\n", iterate->num_close/(double)nsteps);
                    }
                }
            }
            else
            {
                gmx_fatal(FARGS, "Could not converge NPT constraints\n");
            }
        }
    }

    iterate->xprev = iterate->x;
    iterate->x     = *newf;
    iterate->fprev = iterate->f;
    iterate->iter_i++;

    return FALSE;
}
