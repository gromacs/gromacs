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

#include <math.h>
#include "main.h"
#include "constr.h"
#include "copyrite.h"
#include "physics.h"
#include "vec.h"
#include "pbc.h"
#include "smalloc.h"
#include "mdrun.h"
#include "nrnb.h"
#include "domdec.h"
#include "partdec.h"
#include "mtop_util.h"
#include "gmxfio.h"
#include "gmx_omp_nthreads.h"
#include "gmx_omp.h"

typedef struct {
    int    b0;         /* first constraint for this thread */
    int    b1;         /* b1-1 is the last constraint for this thread */
    int    nind;       /* number of indices */
    int   *ind;        /* constraint index for updating atom data */
    int    nind_r;     /* number of indices */
    int   *ind_r;      /* constraint index for updating atom data */
    int    ind_nalloc; /* allocation size of ind and ind_r */
    tensor vir_r_m_dr; /* temporary variable for virial calculation */
} lincs_thread_t;

typedef struct gmx_lincsdata {
    int             ncg;          /* the global number of constraints */
    int             ncg_flex;     /* the global number of flexible constraints */
    int             ncg_triangle; /* the global number of constraints in triangles */
    int             nIter;        /* the number of iterations */
    int             nOrder;       /* the order of the matrix expansion */
    int             nc;           /* the number of constraints */
    int             nc_alloc;     /* the number we allocated memory for */
    int             ncc;          /* the number of constraint connections */
    int             ncc_alloc;    /* the number we allocated memory for */
    real            matlam;       /* the FE lambda value used for filling blc and blmf */
    real           *bllen0;       /* the reference distance in topology A */
    real           *ddist;        /* the reference distance in top B - the r.d. in top A */
    int            *bla;          /* the atom pairs involved in the constraints */
    real           *blc;          /* 1/sqrt(invmass1 + invmass2) */
    real           *blc1;         /* as blc, but with all masses 1 */
    int            *blnr;         /* index into blbnb and blmf */
    int            *blbnb;        /* list of constraint connections */
    int             ntriangle;    /* the local number of constraints in triangles */
    int            *triangle;     /* the list of triangle constraints */
    int            *tri_bits;     /* the bits tell if the matrix element should be used */
    int             ncc_triangle; /* the number of constraint connections in triangles */
    gmx_bool        bCommIter;    /* communicate before each LINCS interation */
    real           *blmf;         /* matrix of mass factors for constraint connections */
    real           *blmf1;        /* as blmf, but with all masses 1 */
    real           *bllen;        /* the reference bond length */
    int             nth;          /* The number of threads doing LINCS */
    lincs_thread_t *th;           /* LINCS thread division */
    unsigned       *atf;          /* atom flags for thread parallelization */
    int             atf_nalloc;   /* allocation size of atf */
    /* arrays for temporary storage in the LINCS algorithm */
    rvec           *tmpv;
    real           *tmpncc;
    real           *tmp1;
    real           *tmp2;
    real           *tmp3;
    real           *tmp4;
    real           *mlambda; /* the Lagrange multipliers * -1 */
    /* storage for the constraint RMS relative deviation output */
    real            rmsd_data[3];
} t_gmx_lincsdata;

real *lincs_rmsd_data(struct gmx_lincsdata *lincsd)
{
    return lincsd->rmsd_data;
}

real lincs_rmsd(struct gmx_lincsdata *lincsd, gmx_bool bSD2)
{
    if (lincsd->rmsd_data[0] > 0)
    {
        return sqrt(lincsd->rmsd_data[bSD2 ? 2 : 1]/lincsd->rmsd_data[0]);
    }
    else
    {
        return 0;
    }
}

/* Do a set of nrec LINCS matrix multiplications.
 * This function will return with up to date thread-local
 * constraint data, without an OpenMP barrier.
 */
static void lincs_matrix_expand(const struct gmx_lincsdata *lincsd,
                                int b0, int b1,
                                const real *blcc,
                                real *rhs1, real *rhs2, real *sol)
{
    int        nrec, rec, b, j, n, nr0, nr1;
    real       mvb, *swap;
    int        ntriangle, tb, bits;
    const int *blnr     = lincsd->blnr, *blbnb = lincsd->blbnb;
    const int *triangle = lincsd->triangle, *tri_bits = lincsd->tri_bits;

    ntriangle = lincsd->ntriangle;
    nrec      = lincsd->nOrder;

    for (rec = 0; rec < nrec; rec++)
    {
#pragma omp barrier
        for (b = b0; b < b1; b++)
        {
            mvb = 0;
            for (n = blnr[b]; n < blnr[b+1]; n++)
            {
                j   = blbnb[n];
                mvb = mvb + blcc[n]*rhs1[j];
            }
            rhs2[b] = mvb;
            sol[b]  = sol[b] + mvb;
        }
        swap = rhs1;
        rhs1 = rhs2;
        rhs2 = swap;
    } /* nrec*(ncons+2*nrtot) flops */

    if (ntriangle > 0)
    {
        /* Perform an extra nrec recursions for only the constraints
         * involved in rigid triangles.
         * In this way their accuracy should come close to those of the other
         * constraints, since traingles of constraints can produce eigenvalues
         * around 0.7, while the effective eigenvalue for bond constraints
         * is around 0.4 (and 0.7*0.7=0.5).
         */
        /* We need to copy the temporary array, since only the elements
         * for constraints involved in triangles are updated and then
         * the pointers are swapped. This saving copying the whole arrary.
         * We need barrier as other threads might still be reading from rhs2.
         */
#pragma omp barrier
        for (b = b0; b < b1; b++)
        {
            rhs2[b] = rhs1[b];
        }
#pragma omp barrier
#pragma omp master
        {
            for (rec = 0; rec < nrec; rec++)
            {
                for (tb = 0; tb < ntriangle; tb++)
                {
                    b    = triangle[tb];
                    bits = tri_bits[tb];
                    mvb  = 0;
                    nr0  = blnr[b];
                    nr1  = blnr[b+1];
                    for (n = nr0; n < nr1; n++)
                    {
                        if (bits & (1<<(n-nr0)))
                        {
                            j   = blbnb[n];
                            mvb = mvb + blcc[n]*rhs1[j];
                        }
                    }
                    rhs2[b] = mvb;
                    sol[b]  = sol[b] + mvb;
                }
                swap = rhs1;
                rhs1 = rhs2;
                rhs2 = swap;
            }
        } /* flops count is missing here */

        /* We need a barrier here as the calling routine will continue
         * to operate on the thread-local constraints without barrier.
         */
#pragma omp barrier
    }
}

static void lincs_update_atoms_noind(int ncons, const int *bla,
                                     real prefac,
                                     const real *fac, rvec *r,
                                     const real *invmass,
                                     rvec *x)
{
    int  b, i, j;
    real mvb, im1, im2, tmp0, tmp1, tmp2;

    if (invmass != NULL)
    {
        for (b = 0; b < ncons; b++)
        {
            i        = bla[2*b];
            j        = bla[2*b+1];
            mvb      = prefac*fac[b];
            im1      = invmass[i];
            im2      = invmass[j];
            tmp0     = r[b][0]*mvb;
            tmp1     = r[b][1]*mvb;
            tmp2     = r[b][2]*mvb;
            x[i][0] -= tmp0*im1;
            x[i][1] -= tmp1*im1;
            x[i][2] -= tmp2*im1;
            x[j][0] += tmp0*im2;
            x[j][1] += tmp1*im2;
            x[j][2] += tmp2*im2;
        } /* 16 ncons flops */
    }
    else
    {
        for (b = 0; b < ncons; b++)
        {
            i        = bla[2*b];
            j        = bla[2*b+1];
            mvb      = prefac*fac[b];
            tmp0     = r[b][0]*mvb;
            tmp1     = r[b][1]*mvb;
            tmp2     = r[b][2]*mvb;
            x[i][0] -= tmp0;
            x[i][1] -= tmp1;
            x[i][2] -= tmp2;
            x[j][0] += tmp0;
            x[j][1] += tmp1;
            x[j][2] += tmp2;
        }
    }
}

static void lincs_update_atoms_ind(int ncons, const int *ind, const int *bla,
                                   real prefac,
                                   const real *fac, rvec *r,
                                   const real *invmass,
                                   rvec *x)
{
    int  bi, b, i, j;
    real mvb, im1, im2, tmp0, tmp1, tmp2;

    if (invmass != NULL)
    {
        for (bi = 0; bi < ncons; bi++)
        {
            b        = ind[bi];
            i        = bla[2*b];
            j        = bla[2*b+1];
            mvb      = prefac*fac[b];
            im1      = invmass[i];
            im2      = invmass[j];
            tmp0     = r[b][0]*mvb;
            tmp1     = r[b][1]*mvb;
            tmp2     = r[b][2]*mvb;
            x[i][0] -= tmp0*im1;
            x[i][1] -= tmp1*im1;
            x[i][2] -= tmp2*im1;
            x[j][0] += tmp0*im2;
            x[j][1] += tmp1*im2;
            x[j][2] += tmp2*im2;
        } /* 16 ncons flops */
    }
    else
    {
        for (bi = 0; bi < ncons; bi++)
        {
            b        = ind[bi];
            i        = bla[2*b];
            j        = bla[2*b+1];
            mvb      = prefac*fac[b];
            tmp0     = r[b][0]*mvb;
            tmp1     = r[b][1]*mvb;
            tmp2     = r[b][2]*mvb;
            x[i][0] -= tmp0;
            x[i][1] -= tmp1;
            x[i][2] -= tmp2;
            x[j][0] += tmp0;
            x[j][1] += tmp1;
            x[j][2] += tmp2;
        } /* 16 ncons flops */
    }
}

static void lincs_update_atoms(struct gmx_lincsdata *li, int th,
                               real prefac,
                               const real *fac, rvec *r,
                               const real *invmass,
                               rvec *x)
{
    if (li->nth == 1)
    {
        /* Single thread, we simply update for all constraints */
        lincs_update_atoms_noind(li->nc, li->bla, prefac, fac, r, invmass, x);
    }
    else
    {
        /* Update the atom vector components for our thread local
         * constraints that only access our local atom range.
         * This can be done without a barrier.
         */
        lincs_update_atoms_ind(li->th[th].nind, li->th[th].ind,
                               li->bla, prefac, fac, r, invmass, x);

        if (li->th[li->nth].nind > 0)
        {
            /* Update the constraints that operate on atoms
             * in multiple thread atom blocks on the master thread.
             */
#pragma omp barrier
#pragma omp master
            {
                lincs_update_atoms_ind(li->th[li->nth].nind,
                                       li->th[li->nth].ind,
                                       li->bla, prefac, fac, r, invmass, x);
            }
        }
    }
}

/* LINCS projection, works on derivatives of the coordinates */
static void do_lincsp(rvec *x, rvec *f, rvec *fp, t_pbc *pbc,
                      struct gmx_lincsdata *lincsd, int th,
                      real *invmass,
                      int econq, real *dvdlambda,
                      gmx_bool bCalcVir, tensor rmdf)
{
    int      b0, b1, b, i, j, k, n;
    real     tmp0, tmp1, tmp2, im1, im2, mvb, rlen, len, wfac, lam;
    rvec     dx;
    int     *bla, *blnr, *blbnb;
    rvec    *r;
    real    *blc, *blmf, *blcc, *rhs1, *rhs2, *sol;

    b0 = lincsd->th[th].b0;
    b1 = lincsd->th[th].b1;

    bla    = lincsd->bla;
    r      = lincsd->tmpv;
    blnr   = lincsd->blnr;
    blbnb  = lincsd->blbnb;
    if (econq != econqForce)
    {
        /* Use mass-weighted parameters */
        blc  = lincsd->blc;
        blmf = lincsd->blmf;
    }
    else
    {
        /* Use non mass-weighted parameters */
        blc  = lincsd->blc1;
        blmf = lincsd->blmf1;
    }
    blcc   = lincsd->tmpncc;
    rhs1   = lincsd->tmp1;
    rhs2   = lincsd->tmp2;
    sol    = lincsd->tmp3;

    /* Compute normalized i-j vectors */
    if (pbc)
    {
        for (b = b0; b < b1; b++)
        {
            pbc_dx_aiuc(pbc, x[bla[2*b]], x[bla[2*b+1]], dx);
            unitv(dx, r[b]);
        }
    }
    else
    {
        for (b = b0; b < b1; b++)
        {
            rvec_sub(x[bla[2*b]], x[bla[2*b+1]], dx);
            unitv(dx, r[b]);
        } /* 16 ncons flops */
    }

#pragma omp barrier
    for (b = b0; b < b1; b++)
    {
        tmp0 = r[b][0];
        tmp1 = r[b][1];
        tmp2 = r[b][2];
        i    = bla[2*b];
        j    = bla[2*b+1];
        for (n = blnr[b]; n < blnr[b+1]; n++)
        {
            k       = blbnb[n];
            blcc[n] = blmf[n]*(tmp0*r[k][0] + tmp1*r[k][1] + tmp2*r[k][2]);
        } /* 6 nr flops */
        mvb = blc[b]*(tmp0*(f[i][0] - f[j][0]) +
                      tmp1*(f[i][1] - f[j][1]) +
                      tmp2*(f[i][2] - f[j][2]));
        rhs1[b] = mvb;
        sol[b]  = mvb;
        /* 7 flops */
    }
    /* Together: 23*ncons + 6*nrtot flops */

    lincs_matrix_expand(lincsd, b0, b1, blcc, rhs1, rhs2, sol);
    /* nrec*(ncons+2*nrtot) flops */

    if (econq == econqDeriv_FlexCon)
    {
        /* We only want to constraint the flexible constraints,
         * so we mask out the normal ones by setting sol to 0.
         */
        for (b = b0; b < b1; b++)
        {
            if (!(lincsd->bllen0[b] == 0 && lincsd->ddist[b] == 0))
            {
                sol[b] = 0;
            }
        }
    }

    /* We multiply sol by blc, so we can use lincs_update_atoms for OpenMP */
    for (b = b0; b < b1; b++)
    {
        sol[b] *= blc[b];
    }

    /* When constraining forces, we should not use mass weighting,
     * so we pass invmass=NULL, which results in the use of 1 for all atoms.
     */
    lincs_update_atoms(lincsd, th, 1.0, sol, r,
                       (econq != econqForce) ? invmass : NULL, fp);

    if (dvdlambda != NULL)
    {
#pragma omp barrier
        for (b = b0; b < b1; b++)
        {
            *dvdlambda -= sol[b]*lincsd->ddist[b];
        }
        /* 10 ncons flops */
    }

    if (bCalcVir)
    {
        /* Constraint virial,
         * determines sum r_bond x delta f,
         * where delta f is the constraint correction
         * of the quantity that is being constrained.
         */
        for (b = b0; b < b1; b++)
        {
            mvb = lincsd->bllen[b]*sol[b];
            for (i = 0; i < DIM; i++)
            {
                tmp1 = mvb*r[b][i];
                for (j = 0; j < DIM; j++)
                {
                    rmdf[i][j] += tmp1*r[b][j];
                }
            }
        } /* 23 ncons flops */
    }
}

static void do_lincs(rvec *x, rvec *xp, matrix box, t_pbc *pbc,
                     struct gmx_lincsdata *lincsd, int th,
                     real *invmass,
                     t_commrec *cr,
                     gmx_bool bCalcLambda,
                     real wangle, int *warn,
                     real invdt, rvec *v,
                     gmx_bool bCalcVir, tensor vir_r_m_dr)
{
    int      b0, b1, b, i, j, k, n, iter;
    real     tmp0, tmp1, tmp2, im1, im2, mvb, rlen, len, len2, dlen2, wfac;
    rvec     dx;
    int     *bla, *blnr, *blbnb;
    rvec    *r;
    real    *blc, *blmf, *bllen, *blcc, *rhs1, *rhs2, *sol, *blc_sol, *mlambda;
    int     *nlocat;

    b0 = lincsd->th[th].b0;
    b1 = lincsd->th[th].b1;

    bla     = lincsd->bla;
    r       = lincsd->tmpv;
    blnr    = lincsd->blnr;
    blbnb   = lincsd->blbnb;
    blc     = lincsd->blc;
    blmf    = lincsd->blmf;
    bllen   = lincsd->bllen;
    blcc    = lincsd->tmpncc;
    rhs1    = lincsd->tmp1;
    rhs2    = lincsd->tmp2;
    sol     = lincsd->tmp3;
    blc_sol = lincsd->tmp4;
    mlambda = lincsd->mlambda;

    if (DOMAINDECOMP(cr) && cr->dd->constraints)
    {
        nlocat = dd_constraints_nlocalatoms(cr->dd);
    }
    else if (PARTDECOMP(cr))
    {
        nlocat = pd_constraints_nlocalatoms(cr->pd);
    }
    else
    {
        nlocat = NULL;
    }

    if (pbc)
    {
        /* Compute normalized i-j vectors */
        for (b = b0; b < b1; b++)
        {
            pbc_dx_aiuc(pbc, x[bla[2*b]], x[bla[2*b+1]], dx);
            unitv(dx, r[b]);
        }
#pragma omp barrier
        for (b = b0; b < b1; b++)
        {
            for (n = blnr[b]; n < blnr[b+1]; n++)
            {
                blcc[n] = blmf[n]*iprod(r[b], r[blbnb[n]]);
            }
            pbc_dx_aiuc(pbc, xp[bla[2*b]], xp[bla[2*b+1]], dx);
            mvb     = blc[b]*(iprod(r[b], dx) - bllen[b]);
            rhs1[b] = mvb;
            sol[b]  = mvb;
        }
    }
    else
    {
        /* Compute normalized i-j vectors */
        for (b = b0; b < b1; b++)
        {
            i       = bla[2*b];
            j       = bla[2*b+1];
            tmp0    = x[i][0] - x[j][0];
            tmp1    = x[i][1] - x[j][1];
            tmp2    = x[i][2] - x[j][2];
            rlen    = gmx_invsqrt(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2);
            r[b][0] = rlen*tmp0;
            r[b][1] = rlen*tmp1;
            r[b][2] = rlen*tmp2;
        } /* 16 ncons flops */

#pragma omp barrier
        for (b = b0; b < b1; b++)
        {
            tmp0 = r[b][0];
            tmp1 = r[b][1];
            tmp2 = r[b][2];
            len  = bllen[b];
            i    = bla[2*b];
            j    = bla[2*b+1];
            for (n = blnr[b]; n < blnr[b+1]; n++)
            {
                k       = blbnb[n];
                blcc[n] = blmf[n]*(tmp0*r[k][0] + tmp1*r[k][1] + tmp2*r[k][2]);
            } /* 6 nr flops */
            mvb = blc[b]*(tmp0*(xp[i][0] - xp[j][0]) +
                          tmp1*(xp[i][1] - xp[j][1]) +
                          tmp2*(xp[i][2] - xp[j][2]) - len);
            rhs1[b] = mvb;
            sol[b]  = mvb;
            /* 10 flops */
        }
        /* Together: 26*ncons + 6*nrtot flops */
    }

    lincs_matrix_expand(lincsd, b0, b1, blcc, rhs1, rhs2, sol);
    /* nrec*(ncons+2*nrtot) flops */

    for (b = b0; b < b1; b++)
    {
        mlambda[b] = blc[b]*sol[b];
    }

    /* Update the coordinates */
    lincs_update_atoms(lincsd, th, 1.0, mlambda, r, invmass, xp);

    /*
     ********  Correction for centripetal effects  ********
     */

    wfac = cos(DEG2RAD*wangle);
    wfac = wfac*wfac;

    for (iter = 0; iter < lincsd->nIter; iter++)
    {
        if ((lincsd->bCommIter && DOMAINDECOMP(cr) && cr->dd->constraints) ||
            PARTDECOMP(cr))
        {
#pragma omp barrier
#pragma omp master
            {
                /* Communicate the corrected non-local coordinates */
                if (DOMAINDECOMP(cr))
                {
                    dd_move_x_constraints(cr->dd, box, xp, NULL, FALSE);
                }
                else
                {
                    pd_move_x_constraints(cr, xp, NULL);
                }
            }
        }

#pragma omp barrier
        for (b = b0; b < b1; b++)
        {
            len = bllen[b];
            if (pbc)
            {
                pbc_dx_aiuc(pbc, xp[bla[2*b]], xp[bla[2*b+1]], dx);
            }
            else
            {
                rvec_sub(xp[bla[2*b]], xp[bla[2*b+1]], dx);
            }
            len2  = len*len;
            dlen2 = 2*len2 - norm2(dx);
            if (dlen2 < wfac*len2 && (nlocat == NULL || nlocat[b]))
            {
                *warn = b;
            }
            if (dlen2 > 0)
            {
                mvb = blc[b]*(len - dlen2*gmx_invsqrt(dlen2));
            }
            else
            {
                mvb = blc[b]*len;
            }
            rhs1[b] = mvb;
            sol[b]  = mvb;
        } /* 20*ncons flops */

        lincs_matrix_expand(lincsd, b0, b1, blcc, rhs1, rhs2, sol);
        /* nrec*(ncons+2*nrtot) flops */

        for (b = b0; b < b1; b++)
        {
            mvb         = blc[b]*sol[b];
            blc_sol[b]  = mvb;
            mlambda[b] += mvb;
        }

        /* Update the coordinates */
        lincs_update_atoms(lincsd, th, 1.0, blc_sol, r, invmass, xp);
    }
    /* nit*ncons*(37+9*nrec) flops */

    if (v != NULL)
    {
        /* Update the velocities */
        lincs_update_atoms(lincsd, th, invdt, mlambda, r, invmass, v);
        /* 16 ncons flops */
    }

    if (nlocat != NULL && bCalcLambda)
    {
        /* In lincs_update_atoms thread might cross-read mlambda */
#pragma omp barrier

        /* Only account for local atoms */
        for (b = b0; b < b1; b++)
        {
            mlambda[b] *= 0.5*nlocat[b];
        }
    }

    if (bCalcVir)
    {
        /* Constraint virial */
        for (b = b0; b < b1; b++)
        {
            tmp0 = -bllen[b]*mlambda[b];
            for (i = 0; i < DIM; i++)
            {
                tmp1 = tmp0*r[b][i];
                for (j = 0; j < DIM; j++)
                {
                    vir_r_m_dr[i][j] -= tmp1*r[b][j];
                }
            }
        } /* 22 ncons flops */
    }

    /* Total:
     * 26*ncons + 6*nrtot + nrec*(ncons+2*nrtot)
     * + nit * (20*ncons + nrec*(ncons+2*nrtot) + 17 ncons)
     *
     * (26+nrec)*ncons + (6+2*nrec)*nrtot
     * + nit * ((37+nrec)*ncons + 2*nrec*nrtot)
     * if nit=1
     * (63+nrec)*ncons + (6+4*nrec)*nrtot
     */
}

void set_lincs_matrix(struct gmx_lincsdata *li, real *invmass, real lambda)
{
    int        i, a1, a2, n, k, sign, center;
    int        end, nk, kk;
    const real invsqrt2 = 0.7071067811865475244;

    for (i = 0; (i < li->nc); i++)
    {
        a1          = li->bla[2*i];
        a2          = li->bla[2*i+1];
        li->blc[i]  = gmx_invsqrt(invmass[a1] + invmass[a2]);
        li->blc1[i] = invsqrt2;
    }

    /* Construct the coupling coefficient matrix blmf */
    li->ntriangle    = 0;
    li->ncc_triangle = 0;
    for (i = 0; (i < li->nc); i++)
    {
        a1 = li->bla[2*i];
        a2 = li->bla[2*i+1];
        for (n = li->blnr[i]; (n < li->blnr[i+1]); n++)
        {
            k = li->blbnb[n];
            if (a1 == li->bla[2*k] || a2 == li->bla[2*k+1])
            {
                sign = -1;
            }
            else
            {
                sign = 1;
            }
            if (a1 == li->bla[2*k] || a1 == li->bla[2*k+1])
            {
                center = a1;
                end    = a2;
            }
            else
            {
                center = a2;
                end    = a1;
            }
            li->blmf[n]  = sign*invmass[center]*li->blc[i]*li->blc[k];
            li->blmf1[n] = sign*0.5;
            if (li->ncg_triangle > 0)
            {
                /* Look for constraint triangles */
                for (nk = li->blnr[k]; (nk < li->blnr[k+1]); nk++)
                {
                    kk = li->blbnb[nk];
                    if (kk != i && kk != k &&
                        (li->bla[2*kk] == end || li->bla[2*kk+1] == end))
                    {
                        if (li->ntriangle == 0 ||
                            li->triangle[li->ntriangle-1] < i)
                        {
                            /* Add this constraint to the triangle list */
                            li->triangle[li->ntriangle] = i;
                            li->tri_bits[li->ntriangle] = 0;
                            li->ntriangle++;
                            if (li->blnr[i+1] - li->blnr[i] > sizeof(li->tri_bits[0])*8 - 1)
                            {
                                gmx_fatal(FARGS, "A constraint is connected to %d constraints, this is more than the %d allowed for constraints participating in triangles",
                                          li->blnr[i+1] - li->blnr[i],
                                          sizeof(li->tri_bits[0])*8-1);
                            }
                        }
                        li->tri_bits[li->ntriangle-1] |= (1<<(n-li->blnr[i]));
                        li->ncc_triangle++;
                    }
                }
            }
        }
    }

    if (debug)
    {
        fprintf(debug, "Of the %d constraints %d participate in triangles\n",
                li->nc, li->ntriangle);
        fprintf(debug, "There are %d couplings of which %d in triangles\n",
                li->ncc, li->ncc_triangle);
    }

    /* Set matlam,
     * so we know with which lambda value the masses have been set.
     */
    li->matlam = lambda;
}

static int count_triangle_constraints(t_ilist *ilist, t_blocka *at2con)
{
    int      ncon1, ncon_tot;
    int      c0, a00, a01, n1, c1, a10, a11, ac1, n2, c2, a20, a21;
    int      ncon_triangle;
    gmx_bool bTriangle;
    t_iatom *ia1, *ia2, *iap;

    ncon1    = ilist[F_CONSTR].nr/3;
    ncon_tot = ncon1 + ilist[F_CONSTRNC].nr/3;

    ia1 = ilist[F_CONSTR].iatoms;
    ia2 = ilist[F_CONSTRNC].iatoms;

    ncon_triangle = 0;
    for (c0 = 0; c0 < ncon_tot; c0++)
    {
        bTriangle = FALSE;
        iap       = constr_iatomptr(ncon1, ia1, ia2, c0);
        a00       = iap[1];
        a01       = iap[2];
        for (n1 = at2con->index[a01]; n1 < at2con->index[a01+1]; n1++)
        {
            c1 = at2con->a[n1];
            if (c1 != c0)
            {
                iap = constr_iatomptr(ncon1, ia1, ia2, c1);
                a10 = iap[1];
                a11 = iap[2];
                if (a10 == a01)
                {
                    ac1 = a11;
                }
                else
                {
                    ac1 = a10;
                }
                for (n2 = at2con->index[ac1]; n2 < at2con->index[ac1+1]; n2++)
                {
                    c2 = at2con->a[n2];
                    if (c2 != c0 && c2 != c1)
                    {
                        iap = constr_iatomptr(ncon1, ia1, ia2, c2);
                        a20 = iap[1];
                        a21 = iap[2];
                        if (a20 == a00 || a21 == a00)
                        {
                            bTriangle = TRUE;
                        }
                    }
                }
            }
        }
        if (bTriangle)
        {
            ncon_triangle++;
        }
    }

    return ncon_triangle;
}

static gmx_bool more_than_two_sequential_constraints(const t_ilist  *ilist,
                                                     const t_blocka *at2con)
{
    t_iatom  *ia1, *ia2, *iap;
    int       ncon1, ncon_tot, c;
    int       a1, a2;
    gmx_bool  bMoreThanTwoSequentialConstraints;

    ncon1    = ilist[F_CONSTR].nr/3;
    ncon_tot = ncon1 + ilist[F_CONSTRNC].nr/3;

    ia1 = ilist[F_CONSTR].iatoms;
    ia2 = ilist[F_CONSTRNC].iatoms;

    bMoreThanTwoSequentialConstraints = FALSE;
    for (c = 0; c < ncon_tot && !bMoreThanTwoSequentialConstraints; c++)
    {
        iap = constr_iatomptr(ncon1, ia1, ia2, c);
        a1  = iap[1];
        a2  = iap[2];
        /* Check if this constraint has constraints connected at both atoms */
        if (at2con->index[a1+1] - at2con->index[a1] > 1 &&
            at2con->index[a2+1] - at2con->index[a2] > 1)
        {
            bMoreThanTwoSequentialConstraints = TRUE;
        }
    }

    return bMoreThanTwoSequentialConstraints;
}

static int int_comp(const void *a, const void *b)
{
    return (*(int *)a) - (*(int *)b);
}

gmx_lincsdata_t init_lincs(FILE *fplog, gmx_mtop_t *mtop,
                           int nflexcon_global, t_blocka *at2con,
                           gmx_bool bPLINCS, int nIter, int nProjOrder)
{
    struct gmx_lincsdata *li;
    int                   mb;
    gmx_moltype_t        *molt;

    if (fplog)
    {
        fprintf(fplog, "\nInitializing%s LINear Constraint Solver\n",
                bPLINCS ? " Parallel" : "");
    }

    snew(li, 1);

    li->ncg      =
        gmx_mtop_ftype_count(mtop, F_CONSTR) +
        gmx_mtop_ftype_count(mtop, F_CONSTRNC);
    li->ncg_flex = nflexcon_global;

    li->nIter  = nIter;
    li->nOrder = nProjOrder;

    li->ncg_triangle = 0;
    li->bCommIter    = FALSE;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molt              = &mtop->moltype[mtop->molblock[mb].type];
        li->ncg_triangle +=
            mtop->molblock[mb].nmol*
            count_triangle_constraints(molt->ilist,
                                       &at2con[mtop->molblock[mb].type]);
        if (bPLINCS && li->bCommIter == FALSE)
        {
            /* Check if we need to communicate not only before LINCS,
             * but also before each iteration.
             * The check for only two sequential constraints is only
             * useful for the common case of H-bond only constraints.
             * With more effort we could also make it useful for small
             * molecules with nr. sequential constraints <= nOrder-1.
             */
            li->bCommIter = (li->nOrder < 1 || more_than_two_sequential_constraints(molt->ilist, &at2con[mtop->molblock[mb].type]));
        }
    }
    if (debug && bPLINCS)
    {
        fprintf(debug, "PLINCS communication before each iteration: %d\n",
                li->bCommIter);
    }

    /* LINCS can run on any number of threads.
     * Currently the number is fixed for the whole simulation,
     * but it could be set in set_lincs().
     */
    li->nth = gmx_omp_nthreads_get(emntLINCS);
    if (li->nth == 1)
    {
        snew(li->th, 1);
    }
    else
    {
        /* Allocate an extra elements for "thread-overlap" constraints */
        snew(li->th, li->nth+1);
    }
    if (debug)
    {
        fprintf(debug, "LINCS: using %d threads\n", li->nth);
    }

    if (bPLINCS || li->ncg_triangle > 0)
    {
        please_cite(fplog, "Hess2008a");
    }
    else
    {
        please_cite(fplog, "Hess97a");
    }

    if (fplog)
    {
        fprintf(fplog, "The number of constraints is %d\n", li->ncg);
        if (bPLINCS)
        {
            fprintf(fplog, "There are inter charge-group constraints,\n"
                    "will communicate selected coordinates each lincs iteration\n");
        }
        if (li->ncg_triangle > 0)
        {
            fprintf(fplog,
                    "%d constraints are involved in constraint triangles,\n"
                    "will apply an additional matrix expansion of order %d for couplings\n"
                    "between constraints inside triangles\n",
                    li->ncg_triangle, li->nOrder);
        }
    }

    return li;
}

/* Sets up the work division over the threads */
static void lincs_thread_setup(struct gmx_lincsdata *li, int natoms)
{
    lincs_thread_t *li_m;
    int             th;
    unsigned       *atf;
    int             a;

    if (natoms > li->atf_nalloc)
    {
        li->atf_nalloc = over_alloc_large(natoms);
        srenew(li->atf, li->atf_nalloc);
    }

    atf = li->atf;
    /* Clear the atom flags */
    for (a = 0; a < natoms; a++)
    {
        atf[a] = 0;
    }

    for (th = 0; th < li->nth; th++)
    {
        lincs_thread_t *li_th;
        int             b;

        li_th = &li->th[th];

        /* The constraints are divided equally over the threads */
        li_th->b0 = (li->nc* th   )/li->nth;
        li_th->b1 = (li->nc*(th+1))/li->nth;

        if (th < sizeof(*atf)*8)
        {
            /* For each atom set a flag for constraints from each */
            for (b = li_th->b0; b < li_th->b1; b++)
            {
                atf[li->bla[b*2]  ] |= (1U<<th);
                atf[li->bla[b*2+1]] |= (1U<<th);
            }
        }
    }

#pragma omp parallel for num_threads(li->nth) schedule(static)
    for (th = 0; th < li->nth; th++)
    {
        lincs_thread_t *li_th;
        unsigned        mask;
        int             b;

        li_th = &li->th[th];

        if (li_th->b1 - li_th->b0 > li_th->ind_nalloc)
        {
            li_th->ind_nalloc = over_alloc_large(li_th->b1-li_th->b0);
            srenew(li_th->ind, li_th->ind_nalloc);
            srenew(li_th->ind_r, li_th->ind_nalloc);
        }

        if (th < sizeof(*atf)*8)
        {
            mask = (1U<<th) - 1U;

            li_th->nind   = 0;
            li_th->nind_r = 0;
            for (b = li_th->b0; b < li_th->b1; b++)
            {
                /* We let the constraint with the lowest thread index
                 * operate on atoms with constraints from multiple threads.
                 */
                if (((atf[li->bla[b*2]]   & mask) == 0) &&
                    ((atf[li->bla[b*2+1]] & mask) == 0))
                {
                    /* Add the constraint to the local atom update index */
                    li_th->ind[li_th->nind++] = b;
                }
                else
                {
                    /* Add the constraint to the rest block */
                    li_th->ind_r[li_th->nind_r++] = b;
                }
            }
        }
        else
        {
            /* We are out of bits, assign all constraints to rest */
            for (b = li_th->b0; b < li_th->b1; b++)
            {
                li_th->ind_r[li_th->nind_r++] = b;
            }
        }
    }

    /* We need to copy all constraints which have not be assigned
     * to a thread to a separate list which will be handled by one thread.
     */
    li_m = &li->th[li->nth];

    li_m->nind = 0;
    for (th = 0; th < li->nth; th++)
    {
        lincs_thread_t *li_th;
        int             b;

        li_th   = &li->th[th];

        if (li_m->nind + li_th->nind_r > li_m->ind_nalloc)
        {
            li_m->ind_nalloc = over_alloc_large(li_m->nind+li_th->nind_r);
            srenew(li_m->ind, li_m->ind_nalloc);
        }

        for (b = 0; b < li_th->nind_r; b++)
        {
            li_m->ind[li_m->nind++] = li_th->ind_r[b];
        }

        if (debug)
        {
            fprintf(debug, "LINCS thread %d: %d constraints\n",
                    th, li_th->nind);
        }
    }

    if (debug)
    {
        fprintf(debug, "LINCS thread r: %d constraints\n",
                li_m->nind);
    }
}


void set_lincs(t_idef *idef, t_mdatoms *md,
               gmx_bool bDynamics, t_commrec *cr,
               struct gmx_lincsdata *li)
{
    int          start, natoms, nflexcon;
    t_blocka     at2con;
    t_iatom     *iatom;
    int          i, k, ncc_alloc, ni, con, nconnect, concon;
    int          type, a1, a2;
    real         lenA = 0, lenB;
    gmx_bool     bLocal;

    li->nc  = 0;
    li->ncc = 0;
    /* Zero the thread index ranges.
     * Otherwise without local constraints we could return with old ranges.
     */
    for (i = 0; i < li->nth; i++)
    {
        li->th[i].b0   = 0;
        li->th[i].b1   = 0;
        li->th[i].nind = 0;
    }
    if (li->nth > 1)
    {
        li->th[li->nth].nind = 0;
    }

    /* This is the local topology, so there are only F_CONSTR constraints */
    if (idef->il[F_CONSTR].nr == 0)
    {
        /* There are no constraints,
         * we do not need to fill any data structures.
         */
        return;
    }

    if (debug)
    {
        fprintf(debug, "Building the LINCS connectivity\n");
    }

    if (DOMAINDECOMP(cr))
    {
        if (cr->dd->constraints)
        {
            dd_get_constraint_range(cr->dd, &start, &natoms);
        }
        else
        {
            natoms = cr->dd->nat_home;
        }
        start = 0;
    }
    else if (PARTDECOMP(cr))
    {
        pd_get_constraint_range(cr->pd, &start, &natoms);
    }
    else
    {
        start  = md->start;
        natoms = md->homenr;
    }
    at2con = make_at2con(start, natoms, idef->il, idef->iparams, bDynamics,
                         &nflexcon);


    if (idef->il[F_CONSTR].nr/3 > li->nc_alloc || li->nc_alloc == 0)
    {
        li->nc_alloc = over_alloc_dd(idef->il[F_CONSTR].nr/3);
        srenew(li->bllen0, li->nc_alloc);
        srenew(li->ddist, li->nc_alloc);
        srenew(li->bla, 2*li->nc_alloc);
        srenew(li->blc, li->nc_alloc);
        srenew(li->blc1, li->nc_alloc);
        srenew(li->blnr, li->nc_alloc+1);
        srenew(li->bllen, li->nc_alloc);
        srenew(li->tmpv, li->nc_alloc);
        srenew(li->tmp1, li->nc_alloc);
        srenew(li->tmp2, li->nc_alloc);
        srenew(li->tmp3, li->nc_alloc);
        srenew(li->tmp4, li->nc_alloc);
        srenew(li->mlambda, li->nc_alloc);
        if (li->ncg_triangle > 0)
        {
            /* This is allocating too much, but it is difficult to improve */
            srenew(li->triangle, li->nc_alloc);
            srenew(li->tri_bits, li->nc_alloc);
        }
    }

    iatom = idef->il[F_CONSTR].iatoms;

    ncc_alloc   = li->ncc_alloc;
    li->blnr[0] = 0;

    ni = idef->il[F_CONSTR].nr/3;

    con           = 0;
    nconnect      = 0;
    li->blnr[con] = nconnect;
    for (i = 0; i < ni; i++)
    {
        bLocal = TRUE;
        type   = iatom[3*i];
        a1     = iatom[3*i+1];
        a2     = iatom[3*i+2];
        lenA   = idef->iparams[type].constr.dA;
        lenB   = idef->iparams[type].constr.dB;
        /* Skip the flexible constraints when not doing dynamics */
        if (bDynamics || lenA != 0 || lenB != 0)
        {
            li->bllen0[con]  = lenA;
            li->ddist[con]   = lenB - lenA;
            /* Set the length to the topology A length */
            li->bllen[con]   = li->bllen0[con];
            li->bla[2*con]   = a1;
            li->bla[2*con+1] = a2;
            /* Construct the constraint connection matrix blbnb */
            for (k = at2con.index[a1-start]; k < at2con.index[a1-start+1]; k++)
            {
                concon = at2con.a[k];
                if (concon != i)
                {
                    if (nconnect >= ncc_alloc)
                    {
                        ncc_alloc = over_alloc_small(nconnect+1);
                        srenew(li->blbnb, ncc_alloc);
                    }
                    li->blbnb[nconnect++] = concon;
                }
            }
            for (k = at2con.index[a2-start]; k < at2con.index[a2-start+1]; k++)
            {
                concon = at2con.a[k];
                if (concon != i)
                {
                    if (nconnect+1 > ncc_alloc)
                    {
                        ncc_alloc = over_alloc_small(nconnect+1);
                        srenew(li->blbnb, ncc_alloc);
                    }
                    li->blbnb[nconnect++] = concon;
                }
            }
            li->blnr[con+1] = nconnect;

            if (cr->dd == NULL)
            {
                /* Order the blbnb matrix to optimize memory access */
                qsort(&(li->blbnb[li->blnr[con]]), li->blnr[con+1]-li->blnr[con],
                      sizeof(li->blbnb[0]), int_comp);
            }
            /* Increase the constraint count */
            con++;
        }
    }

    done_blocka(&at2con);

    /* This is the real number of constraints,
     * without dynamics the flexible constraints are not present.
     */
    li->nc = con;

    li->ncc = li->blnr[con];
    if (cr->dd == NULL)
    {
        /* Since the matrix is static, we can free some memory */
        ncc_alloc = li->ncc;
        srenew(li->blbnb, ncc_alloc);
    }

    if (ncc_alloc > li->ncc_alloc)
    {
        li->ncc_alloc = ncc_alloc;
        srenew(li->blmf, li->ncc_alloc);
        srenew(li->blmf1, li->ncc_alloc);
        srenew(li->tmpncc, li->ncc_alloc);
    }

    if (debug)
    {
        fprintf(debug, "Number of constraints is %d, couplings %d\n",
                li->nc, li->ncc);
    }

    if (li->nth == 1)
    {
        li->th[0].b0 = 0;
        li->th[0].b1 = li->nc;
    }
    else
    {
        lincs_thread_setup(li, md->nr);
    }

    set_lincs_matrix(li, md->invmass, md->lambda);
}

static void lincs_warning(FILE *fplog,
                          gmx_domdec_t *dd, rvec *x, rvec *xprime, t_pbc *pbc,
                          int ncons, int *bla, real *bllen, real wangle,
                          int maxwarn, int *warncount)
{
    int  b, i, j;
    rvec v0, v1;
    real wfac, d0, d1, cosine;
    char buf[STRLEN];

    wfac = cos(DEG2RAD*wangle);

    sprintf(buf, "bonds that rotated more than %g degrees:\n"
            " atom 1 atom 2  angle  previous, current, constraint length\n",
            wangle);
    fprintf(stderr, "%s", buf);
    if (fplog)
    {
        fprintf(fplog, "%s", buf);
    }

    for (b = 0; b < ncons; b++)
    {
        i = bla[2*b];
        j = bla[2*b+1];
        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[i], x[j], v0);
            pbc_dx_aiuc(pbc, xprime[i], xprime[j], v1);
        }
        else
        {
            rvec_sub(x[i], x[j], v0);
            rvec_sub(xprime[i], xprime[j], v1);
        }
        d0     = norm(v0);
        d1     = norm(v1);
        cosine = iprod(v0, v1)/(d0*d1);
        if (cosine < wfac)
        {
            sprintf(buf, " %6d %6d  %5.1f  %8.4f %8.4f    %8.4f\n",
                    ddglatnr(dd, i), ddglatnr(dd, j),
                    RAD2DEG*acos(cosine), d0, d1, bllen[b]);
            fprintf(stderr, "%s", buf);
            if (fplog)
            {
                fprintf(fplog, "%s", buf);
            }
            if (!gmx_isfinite(d1))
            {
                gmx_fatal(FARGS, "Bond length not finite.");
            }

            (*warncount)++;
        }
    }
    if (*warncount > maxwarn)
    {
        too_many_constraint_warnings(econtLINCS, *warncount);
    }
}

static void cconerr(gmx_domdec_t *dd,
                    int ncons, int *bla, real *bllen, rvec *x, t_pbc *pbc,
                    real *ncons_loc, real *ssd, real *max, int *imax)
{
    real       len, d, ma, ssd2, r2;
    int       *nlocat, count, b, im;
    rvec       dx;

    if (dd && dd->constraints)
    {
        nlocat = dd_constraints_nlocalatoms(dd);
    }
    else
    {
        nlocat = 0;
    }

    ma    = 0;
    ssd2  = 0;
    im    = 0;
    count = 0;
    for (b = 0; b < ncons; b++)
    {
        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[bla[2*b]], x[bla[2*b+1]], dx);
        }
        else
        {
            rvec_sub(x[bla[2*b]], x[bla[2*b+1]], dx);
        }
        r2  = norm2(dx);
        len = r2*gmx_invsqrt(r2);
        d   = fabs(len/bllen[b]-1);
        if (d > ma && (nlocat == NULL || nlocat[b]))
        {
            ma = d;
            im = b;
        }
        if (nlocat == NULL)
        {
            ssd2 += d*d;
            count++;
        }
        else
        {
            ssd2  += nlocat[b]*d*d;
            count += nlocat[b];
        }
    }

    *ncons_loc = (nlocat ? 0.5 : 1)*count;
    *ssd       = (nlocat ? 0.5 : 1)*ssd2;
    *max       = ma;
    *imax      = im;
}

static void dump_conf(gmx_domdec_t *dd, struct gmx_lincsdata *li,
                      t_blocka *at2con,
                      char *name, gmx_bool bAll, rvec *x, matrix box)
{
    char  str[STRLEN];
    FILE *fp;
    int   ac0, ac1, i;

    dd_get_constraint_range(dd, &ac0, &ac1);

    sprintf(str, "%s_%d_%d_%d.pdb", name, dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);
    fp = gmx_fio_fopen(str, "w");
    fprintf(fp, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
            10*norm(box[XX]), 10*norm(box[YY]), 10*norm(box[ZZ]),
            90.0, 90.0, 90.0);
    for (i = 0; i < ac1; i++)
    {
        if (i < dd->nat_home || (bAll && i >= ac0 && i < ac1))
        {
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "ATOM", ddglatnr(dd, i), "C", "ALA", ' ', i+1,
                    10*x[i][XX], 10*x[i][YY], 10*x[i][ZZ],
                    1.0, i < dd->nat_tot ? 0.0 : 1.0);
        }
    }
    if (bAll)
    {
        for (i = 0; i < li->nc; i++)
        {
            fprintf(fp, "CONECT%5d%5d\n",
                    ddglatnr(dd, li->bla[2*i]),
                    ddglatnr(dd, li->bla[2*i+1]));
        }
    }
    gmx_fio_fclose(fp);
}

gmx_bool constrain_lincs(FILE *fplog, gmx_bool bLog, gmx_bool bEner,
                         t_inputrec *ir,
                         gmx_large_int_t step,
                         struct gmx_lincsdata *lincsd, t_mdatoms *md,
                         t_commrec *cr,
                         rvec *x, rvec *xprime, rvec *min_proj,
                         matrix box, t_pbc *pbc,
                         real lambda, real *dvdlambda,
                         real invdt, rvec *v,
                         gmx_bool bCalcVir, tensor vir_r_m_dr,
                         int econq,
                         t_nrnb *nrnb,
                         int maxwarn, int *warncount)
{
    char      buf[STRLEN], buf2[22], buf3[STRLEN];
    int       i, warn, p_imax, error;
    real      ncons_loc, p_ssd, p_max = 0;
    rvec      dx;
    gmx_bool  bOK;

    bOK = TRUE;

    if (lincsd->nc == 0 && cr->dd == NULL)
    {
        if (bLog || bEner)
        {
            lincsd->rmsd_data[0] = 0;
            if (ir->eI == eiSD2 && v == NULL)
            {
                i = 2;
            }
            else
            {
                i = 1;
            }
            lincsd->rmsd_data[i] = 0;
        }

        return bOK;
    }

    if (econq == econqCoord)
    {
        if (ir->efep != efepNO)
        {
            if (md->nMassPerturbed && lincsd->matlam != md->lambda)
            {
                set_lincs_matrix(lincsd, md->invmass, md->lambda);
            }

            for (i = 0; i < lincsd->nc; i++)
            {
                lincsd->bllen[i] = lincsd->bllen0[i] + lambda*lincsd->ddist[i];
            }
        }

        if (lincsd->ncg_flex)
        {
            /* Set the flexible constraint lengths to the old lengths */
            if (pbc != NULL)
            {
                for (i = 0; i < lincsd->nc; i++)
                {
                    if (lincsd->bllen[i] == 0)
                    {
                        pbc_dx_aiuc(pbc, x[lincsd->bla[2*i]], x[lincsd->bla[2*i+1]], dx);
                        lincsd->bllen[i] = norm(dx);
                    }
                }
            }
            else
            {
                for (i = 0; i < lincsd->nc; i++)
                {
                    if (lincsd->bllen[i] == 0)
                    {
                        lincsd->bllen[i] =
                            sqrt(distance2(x[lincsd->bla[2*i]],
                                           x[lincsd->bla[2*i+1]]));
                    }
                }
            }
        }

        if (bLog && fplog)
        {
            cconerr(cr->dd, lincsd->nc, lincsd->bla, lincsd->bllen, xprime, pbc,
                    &ncons_loc, &p_ssd, &p_max, &p_imax);
        }

        /* This warn var can be updated by multiple threads
         * at the same time. But as we only need to detect
         * if a warning occured or not, this is not an issue.
         */
        warn = -1;

        /* The OpenMP parallel region of constrain_lincs for coords */
#pragma omp parallel num_threads(lincsd->nth)
        {
            int th = gmx_omp_get_thread_num();

            clear_mat(lincsd->th[th].vir_r_m_dr);

            do_lincs(x, xprime, box, pbc, lincsd, th,
                     md->invmass, cr,
                     bCalcVir || (ir->efep != efepNO),
                     ir->LincsWarnAngle, &warn,
                     invdt, v, bCalcVir,
                     th == 0 ? vir_r_m_dr : lincsd->th[th].vir_r_m_dr);
        }

        if (ir->efep != efepNO)
        {
            real dt_2, dvdl = 0;

            dt_2 = 1.0/(ir->delta_t*ir->delta_t);
            for (i = 0; (i < lincsd->nc); i++)
            {
                dvdl -= lincsd->mlambda[i]*dt_2*lincsd->ddist[i];
            }
            *dvdlambda += dvdl;
        }

        if (bLog && fplog && lincsd->nc > 0)
        {
            fprintf(fplog, "   Rel. Constraint Deviation:  RMS         MAX     between atoms\n");
            fprintf(fplog, "       Before LINCS          %.6f    %.6f %6d %6d\n",
                    sqrt(p_ssd/ncons_loc), p_max,
                    ddglatnr(cr->dd, lincsd->bla[2*p_imax]),
                    ddglatnr(cr->dd, lincsd->bla[2*p_imax+1]));
        }
        if (bLog || bEner)
        {
            cconerr(cr->dd, lincsd->nc, lincsd->bla, lincsd->bllen, xprime, pbc,
                    &ncons_loc, &p_ssd, &p_max, &p_imax);
            /* Check if we are doing the second part of SD */
            if (ir->eI == eiSD2 && v == NULL)
            {
                i = 2;
            }
            else
            {
                i = 1;
            }
            lincsd->rmsd_data[0] = ncons_loc;
            lincsd->rmsd_data[i] = p_ssd;
        }
        else
        {
            lincsd->rmsd_data[0] = 0;
            lincsd->rmsd_data[1] = 0;
            lincsd->rmsd_data[2] = 0;
        }
        if (bLog && fplog && lincsd->nc > 0)
        {
            fprintf(fplog,
                    "        After LINCS          %.6f    %.6f %6d %6d\n\n",
                    sqrt(p_ssd/ncons_loc), p_max,
                    ddglatnr(cr->dd, lincsd->bla[2*p_imax]),
                    ddglatnr(cr->dd, lincsd->bla[2*p_imax+1]));
        }

        if (warn >= 0)
        {
            if (maxwarn >= 0)
            {
                cconerr(cr->dd, lincsd->nc, lincsd->bla, lincsd->bllen, xprime, pbc,
                        &ncons_loc, &p_ssd, &p_max, &p_imax);
                if (MULTISIM(cr))
                {
                    sprintf(buf3, " in simulation %d", cr->ms->sim);
                }
                else
                {
                    buf3[0] = 0;
                }
                sprintf(buf, "\nStep %s, time %g (ps)  LINCS WARNING%s\n"
                        "relative constraint deviation after LINCS:\n"
                        "rms %.6f, max %.6f (between atoms %d and %d)\n",
                        gmx_step_str(step, buf2), ir->init_t+step*ir->delta_t,
                        buf3,
                        sqrt(p_ssd/ncons_loc), p_max,
                        ddglatnr(cr->dd, lincsd->bla[2*p_imax]),
                        ddglatnr(cr->dd, lincsd->bla[2*p_imax+1]));
                if (fplog)
                {
                    fprintf(fplog, "%s", buf);
                }
                fprintf(stderr, "%s", buf);
                lincs_warning(fplog, cr->dd, x, xprime, pbc,
                              lincsd->nc, lincsd->bla, lincsd->bllen,
                              ir->LincsWarnAngle, maxwarn, warncount);
            }
            bOK = (p_max < 0.5);
        }

        if (lincsd->ncg_flex)
        {
            for (i = 0; (i < lincsd->nc); i++)
            {
                if (lincsd->bllen0[i] == 0 && lincsd->ddist[i] == 0)
                {
                    lincsd->bllen[i] = 0;
                }
            }
        }
    }
    else
    {
        /* The OpenMP parallel region of constrain_lincs for derivatives */
#pragma omp parallel num_threads(lincsd->nth)
        {
            int th = gmx_omp_get_thread_num();

            do_lincsp(x, xprime, min_proj, pbc, lincsd, th,
                      md->invmass, econq, ir->efep != efepNO ? dvdlambda : NULL,
                      bCalcVir, th == 0 ? vir_r_m_dr : lincsd->th[th].vir_r_m_dr);
        }
    }

    if (bCalcVir && lincsd->nth > 1)
    {
        for (i = 1; i < lincsd->nth; i++)
        {
            m_add(vir_r_m_dr, lincsd->th[i].vir_r_m_dr, vir_r_m_dr);
        }
    }

    /* count assuming nit=1 */
    inc_nrnb(nrnb, eNR_LINCS, lincsd->nc);
    inc_nrnb(nrnb, eNR_LINCSMAT, (2+lincsd->nOrder)*lincsd->ncc);
    if (lincsd->ntriangle > 0)
    {
        inc_nrnb(nrnb, eNR_LINCSMAT, lincsd->nOrder*lincsd->ncc_triangle);
    }
    if (v)
    {
        inc_nrnb(nrnb, eNR_CONSTR_V, lincsd->nc*2);
    }
    if (bCalcVir)
    {
        inc_nrnb(nrnb, eNR_CONSTR_VIR, lincsd->nc);
    }

    return bOK;
}
