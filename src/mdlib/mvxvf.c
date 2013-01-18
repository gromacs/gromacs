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
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sysstuff.h>
#include <string.h>
#include "typedefs.h"
#include "main.h"
#include "mvdata.h"
#include "network.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "symtab.h"
#include "main.h"
#include "typedefs.h"
#include "vec.h"
#include "tgroup.h"
#include "nrnb.h"
#include "partdec.h"

void move_rvecs(const t_commrec *cr, gmx_bool bForward, gmx_bool bSum,
                int left, int right, rvec vecs[], rvec buf[],
                int shift, t_nrnb *nrnb)
{
    int     i, j, j0 = 137, j1 = 391;
    int     cur, nsum;
    int    *index;
#define next ((cur + 1) % cr->nnodes)
#define prev ((cur - 1 + cr->nnodes) % cr->nnodes)

#define HOMENRI(ind, i) ((ind)[(i)+1] - (ind)[(i)])

    index = pd_index(cr);

    if (bSum)
    {
        cur = (cr->nodeid + pd_shift(cr)) % cr->nnodes;
    }
    else
    {
        cur = cr->nodeid;
    }

    nsum = 0;
    for (i = 0; (i < shift); i++)
    {
        if (bSum)
        {
            if (bForward)
            {
                j0 = index[prev];
                j1 = index[prev+1];
            }
            else
            {
                j0 = index[next];
                j1 = index[next+1];
            }
            for (j = j0; (j < j1); j++)
            {
                clear_rvec(buf[j]);
            }
        }
        /* Forward pulse around the ring, to increasing NODE number */
        if (bForward)
        {
            if (bSum)
            {
                gmx_tx_rx_real(cr,
                               GMX_RIGHT, vecs[index[cur ]], HOMENRI(index, cur )*DIM,
                               GMX_LEFT, buf [index[prev]], HOMENRI(index, prev)*DIM);
            }
            else
            {
                gmx_tx_rx_real(cr,
                               GMX_RIGHT, vecs[index[cur ]], HOMENRI(index, cur )*DIM,
                               GMX_LEFT, vecs[index[prev]], HOMENRI(index, prev)*DIM);
            }
            /* Wait for communication to end */
            gmx_wait(cr, right, left);
        }

        /* Backward pulse around the ring, to decreasing NODE number */
        else
        {
            if (bSum)
            {
                gmx_tx_rx_real(cr,
                               GMX_LEFT, vecs[index[cur ]], HOMENRI(index, cur )*DIM,
                               GMX_RIGHT, buf [index[next]], HOMENRI(index, next)*DIM);
            }
            else
            {
                gmx_tx_rx_real(cr,
                               GMX_LEFT, vecs[index[cur ]], HOMENRI(index, cur )*DIM,
                               GMX_RIGHT, vecs[index[next]], HOMENRI(index, next)*DIM);
            }
            /* Wait for communication to end */
            gmx_wait(cr, left, right);
        }

        /* Actual summation */
        if (bSum)
        {
            for (j = j0; (j < j1); j++)
            {
                rvec_inc(vecs[j], buf[j]);
            }
            nsum += (j1-j0);
        }
        if (bForward)
        {
            cur = prev;
        }
        else
        {
            cur = next;
        }
    }
    if (nsum > 0)
    {
        inc_nrnb(nrnb, eNR_FSUM, nsum);
    }
#undef next
#undef prev
}


void move_reals(const t_commrec *cr, gmx_bool bForward, gmx_bool bSum,
                int left, int right, real reals[], real buf[],
                int shift, t_nrnb *nrnb)
{
    int     i, j, j0 = 137, j1 = 391;
    int     cur, nsum;
    int    *index;
#define next ((cur + 1) % cr->nnodes)
#define prev ((cur - 1 + cr->nnodes) % cr->nnodes)

#define HOMENRI(ind, i) ((ind)[(i)+1] - (ind)[(i)])

    index = pd_index(cr);

    if (bSum)
    {
        cur = (cr->nodeid + pd_shift(cr)) % cr->nnodes;
    }
    else
    {
        cur = cr->nodeid;
    }

    nsum = 0;
    for (i = 0; (i < shift); i++)
    {
        if (bSum)
        {
            if (bForward)
            {
                j0 = index[prev];
                j1 = index[prev+1];
            }
            else
            {
                j0 = index[next];
                j1 = index[next+1];
            }
            for (j = j0; (j < j1); j++)
            {
                buf[j] = 0.0;
            }
        }
        /* Forward pulse around the ring, to increasing NODE number */
        if (bForward)
        {
            if (bSum)
            {
                gmx_tx_rx_real(cr,
                               GMX_RIGHT, reals+index[cur ], HOMENRI(index, cur ),
                               GMX_LEFT, buf+index[prev], HOMENRI(index, prev));
            }
            else
            {
                gmx_tx_rx_real(cr,
                               GMX_RIGHT, reals+index[cur ], HOMENRI(index, cur ),
                               GMX_LEFT, reals+index[prev], HOMENRI(index, prev));
            }
            /* Wait for communication to end */
            gmx_wait(cr, right, left);
        }
        else
        {
            /* Backward pulse around the ring, to decreasing NODE number */
            if (bSum)
            {
                gmx_tx_rx_real(cr,
                               GMX_LEFT, reals+index[cur ], HOMENRI(index, cur ),
                               GMX_RIGHT, buf+index[next], HOMENRI(index, next));
            }
            else
            {
                gmx_tx_rx_real(cr,
                               GMX_LEFT, reals+index[cur ], HOMENRI(index, cur ),
                               GMX_RIGHT, reals+index[next], HOMENRI(index, next));
                /* Wait for communication to end */
            }
            gmx_wait(cr, left, right);
        }

        /* Actual summation */
        if (bSum)
        {
            for (j = j0; (j < j1); j++)
            {
                reals[j] += buf[j];
            }
            nsum += (j1-j0);
        }
        if (bForward)
        {
            cur = prev;
        }
        else
        {
            cur = next;
        }
    }

    if (nsum > 0)
    {
        inc_nrnb(nrnb, eNR_FSUM, nsum/3);
    }
#undef next
#undef prev
}

void move_x(FILE *log, const t_commrec *cr,
            int left, int right, rvec x[],
            t_nrnb *nrnb)
{
    move_rvecs(cr, FALSE, FALSE, left, right, x, NULL, pd_shift(cr), nrnb);
    move_rvecs(cr, TRUE, FALSE, left, right, x, NULL, pd_bshift(cr), nrnb);

    where();
}

void move_rborn(FILE *log, const t_commrec *cr,
                int left, int right, real rborn[],
                t_nrnb *nrnb)
{
    move_reals(cr, FALSE, FALSE, left, right, rborn, NULL, pd_shift(cr), nrnb);
    move_reals(cr, TRUE, FALSE, left, right, rborn, NULL, pd_bshift(cr), nrnb);

    where();
}

void move_f(FILE *log, const t_commrec *cr,
            int left, int right, rvec f[], rvec fadd[],
            t_nrnb *nrnb)
{
    move_rvecs(cr, TRUE, TRUE, left, right, f, fadd, pd_shift(cr), nrnb);
    move_rvecs(cr, FALSE, TRUE, left, right, f, fadd, pd_bshift(cr), nrnb);

    where();
}

void move_gpol(FILE *log, const t_commrec *cr,
               int left, int right, real gpol[], real gpol_add[],
               t_nrnb *nrnb)
{
    move_reals(cr, TRUE, TRUE, left, right, gpol, gpol_add, pd_shift(cr), nrnb);
    move_reals(cr, FALSE, TRUE, left, right, gpol, gpol_add, pd_bshift(cr), nrnb);

    where();
}


void move_cgcm(FILE *log, const t_commrec *cr, rvec cg_cm[])
{
    int  i, start, nr;
    int  cur = cr->nodeid;
    int *cgindex;

#define next ((cur+1) % cr->nnodes)

    cgindex = pd_cgindex(cr);

    for (i = 0; (i < cr->nnodes-1); i++)
    {
        start = cgindex[cur];
        nr    = cgindex[cur+1] - start;

        gmx_tx(cr, GMX_LEFT, cg_cm[start], nr*sizeof(cg_cm[0]));
#ifdef DEBUG
        fprintf(log, "move_cgcm: TX start=%d, nr=%d\n", start, nr);
#endif
        start = cgindex[next];
        nr    = cgindex[next+1] - start;

        gmx_rx(cr, GMX_RIGHT, cg_cm[start], nr*sizeof(cg_cm[0]));
#ifdef DEBUG
        fprintf(log, "move_cgcm: RX start=%d, nr=%d\n", start, nr);
#endif
        gmx_tx_wait(cr, GMX_LEFT);
        gmx_rx_wait(cr, GMX_RIGHT);

        if (debug)
        {
            fprintf(debug, "cgcm[0][XX] %f\n", cg_cm[0][XX]);
        }

        cur = next;
    }
#undef next
}
