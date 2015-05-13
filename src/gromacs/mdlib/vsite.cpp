/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/vsite.h"

#include <stdio.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

/* Routines to send/recieve coordinates and force
 * of constructing atoms.
 */

static int pbc_rvec_sub(const t_pbc *pbc, const rvec xi, const rvec xj, rvec dx)
{
    if (pbc)
    {
        return pbc_dx_aiuc(pbc, xi, xj, dx);
    }
    else
    {
        rvec_sub(xi, xj, dx);
        return CENTRAL;
    }
}

/* Vsite construction routines */

static void constr_vsite2(rvec xi, rvec xj, rvec x, real a, t_pbc *pbc)
{
    real b;
    rvec dx;

    b = 1.0-a;
    /* 1 flop */

    if (pbc)
    {
        pbc_dx_aiuc(pbc, xj, xi, dx);
        x[XX] = xi[XX] + a*dx[XX];
        x[YY] = xi[YY] + a*dx[YY];
        x[ZZ] = xi[ZZ] + a*dx[ZZ];
    }
    else
    {
        x[XX] = b*xi[XX] + a*xj[XX];
        x[YY] = b*xi[YY] + a*xj[YY];
        x[ZZ] = b*xi[ZZ] + a*xj[ZZ];
        /* 9 Flops */
    }

    /* TOTAL: 10 flops */
}

static void constr_vsite3(rvec xi, rvec xj, rvec xk, rvec x, real a, real b,
                          t_pbc *pbc)
{
    real c;
    rvec dxj, dxk;

    c = 1.0-a-b;
    /* 2 flops */

    if (pbc)
    {
        pbc_dx_aiuc(pbc, xj, xi, dxj);
        pbc_dx_aiuc(pbc, xk, xi, dxk);
        x[XX] = xi[XX] + a*dxj[XX] + b*dxk[XX];
        x[YY] = xi[YY] + a*dxj[YY] + b*dxk[YY];
        x[ZZ] = xi[ZZ] + a*dxj[ZZ] + b*dxk[ZZ];
    }
    else
    {
        x[XX] = c*xi[XX] + a*xj[XX] + b*xk[XX];
        x[YY] = c*xi[YY] + a*xj[YY] + b*xk[YY];
        x[ZZ] = c*xi[ZZ] + a*xj[ZZ] + b*xk[ZZ];
        /* 15 Flops */
    }

    /* TOTAL: 17 flops */
}

static void constr_vsite3FD(rvec xi, rvec xj, rvec xk, rvec x, real a, real b,
                            t_pbc *pbc)
{
    rvec xij, xjk, temp;
    real c;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    /* 6 flops */

    /* temp goes from i to a point on the line jk */
    temp[XX] = xij[XX] + a*xjk[XX];
    temp[YY] = xij[YY] + a*xjk[YY];
    temp[ZZ] = xij[ZZ] + a*xjk[ZZ];
    /* 6 flops */

    c = b*gmx_invsqrt(iprod(temp, temp));
    /* 6 + 10 flops */

    x[XX] = xi[XX] + c*temp[XX];
    x[YY] = xi[YY] + c*temp[YY];
    x[ZZ] = xi[ZZ] + c*temp[ZZ];
    /* 6 Flops */

    /* TOTAL: 34 flops */
}

static void constr_vsite3FAD(rvec xi, rvec xj, rvec xk, rvec x, real a, real b, t_pbc *pbc)
{
    rvec xij, xjk, xp;
    real a1, b1, c1, invdij;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    /* 6 flops */

    invdij = gmx_invsqrt(iprod(xij, xij));
    c1     = invdij * invdij * iprod(xij, xjk);
    xp[XX] = xjk[XX] - c1*xij[XX];
    xp[YY] = xjk[YY] - c1*xij[YY];
    xp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
    a1     = a*invdij;
    b1     = b*gmx_invsqrt(iprod(xp, xp));
    /* 45 */

    x[XX] = xi[XX] + a1*xij[XX] + b1*xp[XX];
    x[YY] = xi[YY] + a1*xij[YY] + b1*xp[YY];
    x[ZZ] = xi[ZZ] + a1*xij[ZZ] + b1*xp[ZZ];
    /* 12 Flops */

    /* TOTAL: 63 flops */
}

static void constr_vsite3OUT(rvec xi, rvec xj, rvec xk, rvec x,
                             real a, real b, real c, t_pbc *pbc)
{
    rvec xij, xik, temp;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xi, xik);
    cprod(xij, xik, temp);
    /* 15 Flops */

    x[XX] = xi[XX] + a*xij[XX] + b*xik[XX] + c*temp[XX];
    x[YY] = xi[YY] + a*xij[YY] + b*xik[YY] + c*temp[YY];
    x[ZZ] = xi[ZZ] + a*xij[ZZ] + b*xik[ZZ] + c*temp[ZZ];
    /* 18 Flops */

    /* TOTAL: 33 flops */
}

static void constr_vsite4FD(rvec xi, rvec xj, rvec xk, rvec xl, rvec x,
                            real a, real b, real c, t_pbc *pbc)
{
    rvec xij, xjk, xjl, temp;
    real d;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    pbc_rvec_sub(pbc, xl, xj, xjl);
    /* 9 flops */

    /* temp goes from i to a point on the plane jkl */
    temp[XX] = xij[XX] + a*xjk[XX] + b*xjl[XX];
    temp[YY] = xij[YY] + a*xjk[YY] + b*xjl[YY];
    temp[ZZ] = xij[ZZ] + a*xjk[ZZ] + b*xjl[ZZ];
    /* 12 flops */

    d = c*gmx_invsqrt(iprod(temp, temp));
    /* 6 + 10 flops */

    x[XX] = xi[XX] + d*temp[XX];
    x[YY] = xi[YY] + d*temp[YY];
    x[ZZ] = xi[ZZ] + d*temp[ZZ];
    /* 6 Flops */

    /* TOTAL: 43 flops */
}


static void constr_vsite4FDN(rvec xi, rvec xj, rvec xk, rvec xl, rvec x,
                             real a, real b, real c, t_pbc *pbc)
{
    rvec xij, xik, xil, ra, rb, rja, rjb, rm;
    real d;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xi, xik);
    pbc_rvec_sub(pbc, xl, xi, xil);
    /* 9 flops */

    ra[XX] = a*xik[XX];
    ra[YY] = a*xik[YY];
    ra[ZZ] = a*xik[ZZ];

    rb[XX] = b*xil[XX];
    rb[YY] = b*xil[YY];
    rb[ZZ] = b*xil[ZZ];

    /* 6 flops */

    rvec_sub(ra, xij, rja);
    rvec_sub(rb, xij, rjb);
    /* 6 flops */

    cprod(rja, rjb, rm);
    /* 9 flops */

    d = c*gmx_invsqrt(norm2(rm));
    /* 5+5+1 flops */

    x[XX] = xi[XX] + d*rm[XX];
    x[YY] = xi[YY] + d*rm[YY];
    x[ZZ] = xi[ZZ] + d*rm[ZZ];
    /* 6 Flops */

    /* TOTAL: 47 flops */
}


static int constr_vsiten(t_iatom *ia, t_iparams ip[],
                         rvec *x, t_pbc *pbc)
{
    rvec x1, dx;
    dvec dsum;
    int  n3, av, ai, i;
    real a;

    n3 = 3*ip[ia[0]].vsiten.n;
    av = ia[1];
    ai = ia[2];
    copy_rvec(x[ai], x1);
    clear_dvec(dsum);
    for (i = 3; i < n3; i += 3)
    {
        ai = ia[i+2];
        a  = ip[ia[i]].vsiten.a;
        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[ai], x1, dx);
        }
        else
        {
            rvec_sub(x[ai], x1, dx);
        }
        dsum[XX] += a*dx[XX];
        dsum[YY] += a*dx[YY];
        dsum[ZZ] += a*dx[ZZ];
        /* 9 Flops */
    }

    x[av][XX] = x1[XX] + dsum[XX];
    x[av][YY] = x1[YY] + dsum[YY];
    x[av][ZZ] = x1[ZZ] + dsum[ZZ];

    return n3;
}


void construct_vsites_thread(gmx_vsite_t *vsite,
                             rvec x[],
                             real dt, rvec *v,
                             t_iparams ip[], t_ilist ilist[],
                             t_pbc *pbc_null)
{
    gmx_bool   bPBCAll;
    rvec       xpbc, xv, vv, dx;
    real       a1, b1, c1, inv_dt;
    int        i, inc, nra, nr, tp, ftype;
    t_iatom    avsite, ai, aj, ak, al, pbc_atom;
    t_iatom   *ia;
    t_pbc     *pbc_null2;
    int       *vsite_pbc, ishift;

    if (v != NULL)
    {
        inv_dt = 1.0/dt;
    }
    else
    {
        inv_dt = 1.0;
    }

    bPBCAll = (pbc_null != NULL && !vsite->bHaveChargeGroups);

    pbc_null2 = NULL;
    vsite_pbc = NULL;
    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if ((interaction_function[ftype].flags & IF_VSITE) &&
            ilist[ftype].nr > 0)
        {
            nra    = interaction_function[ftype].nratoms;
            inc    = 1 + nra;
            nr     = ilist[ftype].nr;
            ia     = ilist[ftype].iatoms;

            if (bPBCAll)
            {
                pbc_null2 = pbc_null;
            }
            else if (pbc_null != NULL)
            {
                vsite_pbc = vsite->vsite_pbc_loc[ftype-F_VSITE2];
            }

            for (i = 0; i < nr; )
            {
                tp   = ia[0];

                /* The vsite and constructing atoms */
                avsite = ia[1];
                ai     = ia[2];

                /* Constants for constructing vsites */
                a1   = ip[tp].vsite.a;
                /* Check what kind of pbc we need to use */
                if (bPBCAll)
                {
                    /* No charge groups, vsite follows its own pbc */
                    pbc_atom = avsite;
                    copy_rvec(x[avsite], xpbc);
                }
                else if (vsite_pbc != NULL)
                {
                    pbc_atom = vsite_pbc[i/(1+nra)];
                    if (pbc_atom > -2)
                    {
                        if (pbc_atom >= 0)
                        {
                            /* We need to copy the coordinates here,
                             * single for single atom cg's pbc_atom
                             * is the vsite itself.
                             */
                            copy_rvec(x[pbc_atom], xpbc);
                        }
                        pbc_null2 = pbc_null;
                    }
                    else
                    {
                        pbc_null2 = NULL;
                    }
                }
                else
                {
                    pbc_atom = -2;
                }
                /* Copy the old position */
                copy_rvec(x[avsite], xv);

                /* Construct the vsite depending on type */
                switch (ftype)
                {
                    case F_VSITE2:
                        aj = ia[3];
                        constr_vsite2(x[ai], x[aj], x[avsite], a1, pbc_null2);
                        break;
                    case F_VSITE3:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        constr_vsite3(x[ai], x[aj], x[ak], x[avsite], a1, b1, pbc_null2);
                        break;
                    case F_VSITE3FD:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        constr_vsite3FD(x[ai], x[aj], x[ak], x[avsite], a1, b1, pbc_null2);
                        break;
                    case F_VSITE3FAD:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        constr_vsite3FAD(x[ai], x[aj], x[ak], x[avsite], a1, b1, pbc_null2);
                        break;
                    case F_VSITE3OUT:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        constr_vsite3OUT(x[ai], x[aj], x[ak], x[avsite], a1, b1, c1, pbc_null2);
                        break;
                    case F_VSITE4FD:
                        aj = ia[3];
                        ak = ia[4];
                        al = ia[5];
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        constr_vsite4FD(x[ai], x[aj], x[ak], x[al], x[avsite], a1, b1, c1,
                                        pbc_null2);
                        break;
                    case F_VSITE4FDN:
                        aj = ia[3];
                        ak = ia[4];
                        al = ia[5];
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        constr_vsite4FDN(x[ai], x[aj], x[ak], x[al], x[avsite], a1, b1, c1,
                                         pbc_null2);
                        break;
                    case F_VSITEN:
                        inc = constr_vsiten(ia, ip, x, pbc_null2);
                        break;
                    default:
                        gmx_fatal(FARGS, "No such vsite type %d in %s, line %d",
                                  ftype, __FILE__, __LINE__);
                }

                if (pbc_atom >= 0)
                {
                    /* Match the pbc of this vsite to the rest of its charge group */
                    ishift = pbc_dx_aiuc(pbc_null, x[avsite], xpbc, dx);
                    if (ishift != CENTRAL)
                    {
                        rvec_add(xpbc, dx, x[avsite]);
                    }
                }
                if (v != NULL)
                {
                    /* Calculate velocity of vsite... */
                    rvec_sub(x[avsite], xv, vv);
                    svmul(inv_dt, vv, v[avsite]);
                }

                /* Increment loop variables */
                i  += inc;
                ia += inc;
            }
        }
    }
}

void construct_vsites(gmx_vsite_t *vsite,
                      rvec x[],
                      real dt, rvec *v,
                      t_iparams ip[], t_ilist ilist[],
                      int ePBC, gmx_bool bMolPBC,
                      t_commrec *cr, matrix box)
{
    t_pbc     pbc, *pbc_null;
    gmx_bool  bDomDec;

    bDomDec = cr && DOMAINDECOMP(cr);

    /* We only need to do pbc when we have inter-cg vsites */
    if (ePBC != epbcNONE && (bDomDec || bMolPBC) && vsite->n_intercg_vsite)
    {
        /* This is wasting some CPU time as we now do this multiple times
         * per MD step. But how often do we have vsites with full pbc?
         */
        pbc_null = set_pbc_dd(&pbc, ePBC, cr != NULL ? cr->dd : NULL, FALSE, box);
    }
    else
    {
        pbc_null = NULL;
    }

    if (bDomDec)
    {
        dd_move_x_vsites(cr->dd, box, x);
    }

    if (vsite->nthreads == 1)
    {
        construct_vsites_thread(vsite,
                                x, dt, v,
                                ip, ilist,
                                pbc_null);
    }
    else
    {
#pragma omp parallel num_threads(vsite->nthreads)
        {
            construct_vsites_thread(vsite,
                                    x, dt, v,
                                    ip, vsite->tdata[gmx_omp_get_thread_num()].ilist,
                                    pbc_null);
        }
        /* Now we can construct the vsites that might depend on other vsites */
        construct_vsites_thread(vsite,
                                x, dt, v,
                                ip, vsite->tdata[vsite->nthreads].ilist,
                                pbc_null);
    }
}

static void spread_vsite2(t_iatom ia[], real a,
                          rvec x[], rvec f[], rvec fshift[],
                          t_pbc *pbc, t_graph *g)
{
    rvec    fi, fj, dx;
    t_iatom av, ai, aj;
    ivec    di;
    int     siv, sij;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];

    svmul(1-a, f[av], fi);
    svmul(  a, f[av], fj);
    /* 7 flop */

    rvec_inc(f[ai], fi);
    rvec_inc(f[aj], fj);
    /* 6 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, av), di);
        siv = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), di);
        sij = IVEC2IS(di);
    }
    else if (pbc)
    {
        siv = pbc_dx_aiuc(pbc, x[ai], x[av], dx);
        sij = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
    }
    else
    {
        siv = CENTRAL;
        sij = CENTRAL;
    }

    if (fshift && (siv != CENTRAL || sij != CENTRAL))
    {
        rvec_inc(fshift[siv], f[av]);
        rvec_dec(fshift[CENTRAL], fi);
        rvec_dec(fshift[sij], fj);
    }

    /* TOTAL: 13 flops */
}

void construct_vsites_mtop(gmx_vsite_t *vsite,
                           gmx_mtop_t *mtop, rvec x[])
{
    int             as, mb, mol;
    gmx_molblock_t *molb;
    gmx_moltype_t  *molt;

    as = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];
        for (mol = 0; mol < molb->nmol; mol++)
        {
            construct_vsites(vsite, x+as, 0.0, NULL,
                             mtop->ffparams.iparams, molt->ilist,
                             epbcNONE, TRUE, NULL, NULL);
            as += molt->atoms.nr;
        }
    }
}

static void spread_vsite3(t_iatom ia[], real a, real b,
                          rvec x[], rvec f[], rvec fshift[],
                          t_pbc *pbc, t_graph *g)
{
    rvec    fi, fj, fk, dx;
    atom_id av, ai, aj, ak;
    ivec    di;
    int     siv, sij, sik;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];

    svmul(1-a-b, f[av], fi);
    svmul(    a, f[av], fj);
    svmul(    b, f[av], fk);
    /* 11 flops */

    rvec_inc(f[ai], fi);
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    /* 9 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, ia[1]), di);
        siv = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), di);
        sij = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, ak), di);
        sik = IVEC2IS(di);
    }
    else if (pbc)
    {
        siv = pbc_dx_aiuc(pbc, x[ai], x[av], dx);
        sij = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        sik = pbc_dx_aiuc(pbc, x[ai], x[ak], dx);
    }
    else
    {
        siv = CENTRAL;
        sij = CENTRAL;
        sik = CENTRAL;
    }

    if (fshift && (siv != CENTRAL || sij != CENTRAL || sik != CENTRAL))
    {
        rvec_inc(fshift[siv], f[av]);
        rvec_dec(fshift[CENTRAL], fi);
        rvec_dec(fshift[sij], fj);
        rvec_dec(fshift[sik], fk);
    }

    /* TOTAL: 20 flops */
}

static void spread_vsite3FD(t_iatom ia[], real a, real b,
                            rvec x[], rvec f[], rvec fshift[],
                            gmx_bool VirCorr, matrix dxdf,
                            t_pbc *pbc, t_graph *g)
{
    real    c, invl, fproj, a1;
    rvec    xvi, xij, xjk, xix, fv, temp;
    t_iatom av, ai, aj, ak;
    int     svi, sji, skj;
    ivec    di;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    copy_rvec(f[av], fv);

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    /* 6 flops */

    /* xix goes from i to point x on the line jk */
    xix[XX] = xij[XX]+a*xjk[XX];
    xix[YY] = xij[YY]+a*xjk[YY];
    xix[ZZ] = xij[ZZ]+a*xjk[ZZ];
    /* 6 flops */

    invl = gmx_invsqrt(iprod(xix, xix));
    c    = b*invl;
    /* 4 + ?10? flops */

    fproj = iprod(xix, fv)*invl*invl; /* = (xix . f)/(xix . xix) */

    temp[XX] = c*(fv[XX]-fproj*xix[XX]);
    temp[YY] = c*(fv[YY]-fproj*xix[YY]);
    temp[ZZ] = c*(fv[ZZ]-fproj*xix[ZZ]);
    /* 16 */

    /* c is already calculated in constr_vsite3FD
       storing c somewhere will save 26 flops!     */

    a1         = 1-a;
    f[ai][XX] += fv[XX] - temp[XX];
    f[ai][YY] += fv[YY] - temp[YY];
    f[ai][ZZ] += fv[ZZ] - temp[ZZ];
    f[aj][XX] += a1*temp[XX];
    f[aj][YY] += a1*temp[YY];
    f[aj][ZZ] += a1*temp[ZZ];
    f[ak][XX] += a*temp[XX];
    f[ak][YY] += a*temp[YY];
    f[ak][ZZ] += a*temp[ZZ];
    /* 19 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, aj), di);
        skj = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sji != CENTRAL || skj != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - (1 + a)*temp[XX];
        fshift[CENTRAL][YY] += fv[YY] - (1 + a)*temp[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - (1 + a)*temp[ZZ];
        fshift[    sji][XX] += temp[XX];
        fshift[    sji][YY] += temp[YY];
        fshift[    sji][ZZ] += temp[ZZ];
        fshift[    skj][XX] += a*temp[XX];
        fshift[    skj][YY] += a*temp[YY];
        fshift[    skj][ZZ] += a*temp[ZZ];
    }

    if (VirCorr)
    {
        /* When VirCorr=TRUE, the virial for the current forces is not
         * calculated from the redistributed forces. This means that
         * the effect of non-linear virtual site constructions on the virial
         * needs to be added separately. This contribution can be calculated
         * in many ways, but the simplest and cheapest way is to use
         * the first constructing atom ai as a reference position in space:
         * subtract (xv-xi)*fv and add (xj-xi)*fj + (xk-xi)*fk.
         */
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                /* As xix is a linear combination of j and k, use that here */
                dxdf[i][j] += -xiv[i]*fv[j] + xix[i]*temp[j];
            }
        }
    }

    /* TOTAL: 61 flops */
}

static void spread_vsite3FAD(t_iatom ia[], real a, real b,
                             rvec x[], rvec f[], rvec fshift[],
                             gmx_bool VirCorr, matrix dxdf,
                             t_pbc *pbc, t_graph *g)
{
    rvec    xvi, xij, xjk, xperp, Fpij, Fppp, fv, f1, f2, f3;
    real    a1, b1, c1, c2, invdij, invdij2, invdp, fproj;
    t_iatom av, ai, aj, ak;
    int     svi, sji, skj, d;
    ivec    di;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    copy_rvec(f[ia[1]], fv);

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    /* 6 flops */

    invdij    = gmx_invsqrt(iprod(xij, xij));
    invdij2   = invdij * invdij;
    c1        = iprod(xij, xjk) * invdij2;
    xperp[XX] = xjk[XX] - c1*xij[XX];
    xperp[YY] = xjk[YY] - c1*xij[YY];
    xperp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
    /* xperp in plane ijk, perp. to ij */
    invdp = gmx_invsqrt(iprod(xperp, xperp));
    a1    = a*invdij;
    b1    = b*invdp;
    /* 45 flops */

    /* a1, b1 and c1 are already calculated in constr_vsite3FAD
       storing them somewhere will save 45 flops!     */

    fproj = iprod(xij, fv)*invdij2;
    svmul(fproj,                      xij,  Fpij);    /* proj. f on xij */
    svmul(iprod(xperp, fv)*invdp*invdp, xperp, Fppp); /* proj. f on xperp */
    svmul(b1*fproj,                   xperp, f3);
    /* 23 flops */

    rvec_sub(fv, Fpij, f1); /* f1 = f - Fpij */
    rvec_sub(f1, Fppp, f2); /* f2 = f - Fpij - Fppp */
    for (d = 0; (d < DIM); d++)
    {
        f1[d] *= a1;
        f2[d] *= b1;
    }
    /* 12 flops */

    c2         = 1+c1;
    f[ai][XX] += fv[XX] - f1[XX] + c1*f2[XX] + f3[XX];
    f[ai][YY] += fv[YY] - f1[YY] + c1*f2[YY] + f3[YY];
    f[ai][ZZ] += fv[ZZ] - f1[ZZ] + c1*f2[ZZ] + f3[ZZ];
    f[aj][XX] +=          f1[XX] - c2*f2[XX] - f3[XX];
    f[aj][YY] +=          f1[YY] - c2*f2[YY] - f3[YY];
    f[aj][ZZ] +=          f1[ZZ] - c2*f2[ZZ] - f3[ZZ];
    f[ak][XX] +=                      f2[XX];
    f[ak][YY] +=                      f2[YY];
    f[ak][ZZ] +=                      f2[ZZ];
    /* 30 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, aj), di);
        skj = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sji != CENTRAL || skj != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - f1[XX] - (1-c1)*f2[XX] + f3[XX];
        fshift[CENTRAL][YY] += fv[YY] - f1[YY] - (1-c1)*f2[YY] + f3[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - f1[ZZ] - (1-c1)*f2[ZZ] + f3[ZZ];
        fshift[    sji][XX] +=          f1[XX] -    c1 *f2[XX] - f3[XX];
        fshift[    sji][YY] +=          f1[YY] -    c1 *f2[YY] - f3[YY];
        fshift[    sji][ZZ] +=          f1[ZZ] -    c1 *f2[ZZ] - f3[ZZ];
        fshift[    skj][XX] +=                          f2[XX];
        fshift[    skj][YY] +=                          f2[YY];
        fshift[    skj][ZZ] +=                          f2[ZZ];
    }

    if (VirCorr)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                /* Note that xik=xij+xjk, so we have to add xij*f2 */
                dxdf[i][j] +=
                    -xiv[i]*fv[j]
                    + xij[i]*(f1[j] + (1 - c2)*f2[j] - f3[j])
                    + xjk[i]*f2[j];
            }
        }
    }

    /* TOTAL: 113 flops */
}

static void spread_vsite3OUT(t_iatom ia[], real a, real b, real c,
                             rvec x[], rvec f[], rvec fshift[],
                             gmx_bool VirCorr, matrix dxdf,
                             t_pbc *pbc, t_graph *g)
{
    rvec    xvi, xij, xik, fv, fj, fk;
    real    cfx, cfy, cfz;
    atom_id av, ai, aj, ak;
    ivec    di;
    int     svi, sji, ski;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    ski = pbc_rvec_sub(pbc, x[ak], x[ai], xik);
    /* 6 Flops */

    copy_rvec(f[av], fv);

    cfx = c*fv[XX];
    cfy = c*fv[YY];
    cfz = c*fv[ZZ];
    /* 3 Flops */

    fj[XX] = a*fv[XX]     -  xik[ZZ]*cfy +  xik[YY]*cfz;
    fj[YY] =  xik[ZZ]*cfx + a*fv[YY]     -  xik[XX]*cfz;
    fj[ZZ] = -xik[YY]*cfx +  xik[XX]*cfy + a*fv[ZZ];

    fk[XX] = b*fv[XX]     +  xij[ZZ]*cfy -  xij[YY]*cfz;
    fk[YY] = -xij[ZZ]*cfx + b*fv[YY]     +  xij[XX]*cfz;
    fk[ZZ] =  xij[YY]*cfx -  xij[XX]*cfy + b*fv[ZZ];
    /* 30 Flops */

    f[ai][XX] += fv[XX] - fj[XX] - fk[XX];
    f[ai][YY] += fv[YY] - fj[YY] - fk[YY];
    f[ai][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ];
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    /* 15 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, ai), di);
        ski = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sji != CENTRAL || ski != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - fj[XX] - fk[XX];
        fshift[CENTRAL][YY] += fv[YY] - fj[YY] - fk[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ];
        rvec_inc(fshift[sji], fj);
        rvec_inc(fshift[ski], fk);
    }

    if (VirCorr)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i]*fv[j] + xij[i]*fj[j] + xik[i]*fk[j];
            }
        }
    }

    /* TOTAL: 54 flops */
}

static void spread_vsite4FD(t_iatom ia[], real a, real b, real c,
                            rvec x[], rvec f[], rvec fshift[],
                            gmx_bool VirCorr, matrix dxdf,
                            t_pbc *pbc, t_graph *g)
{
    real    d, invl, fproj, a1;
    rvec    xvi, xij, xjk, xjl, xix, fv, temp;
    atom_id av, ai, aj, ak, al;
    ivec    di;
    int     svi, sji, skj, slj, m;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    al = ia[5];

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    slj = pbc_rvec_sub(pbc, x[al], x[aj], xjl);
    /* 9 flops */

    /* xix goes from i to point x on the plane jkl */
    for (m = 0; m < DIM; m++)
    {
        xix[m] = xij[m] + a*xjk[m] + b*xjl[m];
    }
    /* 12 flops */

    invl = gmx_invsqrt(iprod(xix, xix));
    d    = c*invl;
    /* 4 + ?10? flops */

    copy_rvec(f[av], fv);

    fproj = iprod(xix, fv)*invl*invl; /* = (xix . f)/(xix . xix) */

    for (m = 0; m < DIM; m++)
    {
        temp[m] = d*(fv[m] - fproj*xix[m]);
    }
    /* 16 */

    /* c is already calculated in constr_vsite3FD
       storing c somewhere will save 35 flops!     */

    a1 = 1 - a - b;
    for (m = 0; m < DIM; m++)
    {
        f[ai][m] += fv[m] - temp[m];
        f[aj][m] += a1*temp[m];
        f[ak][m] += a*temp[m];
        f[al][m] += b*temp[m];
    }
    /* 26 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, aj), di);
        skj = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, al), SHIFT_IVEC(g, aj), di);
        slj = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift &&
        (svi != CENTRAL || sji != CENTRAL || skj != CENTRAL || slj != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        for (m = 0; m < DIM; m++)
        {
            fshift[CENTRAL][m] += fv[m] - (1 + a + b)*temp[m];
            fshift[    sji][m] += temp[m];
            fshift[    skj][m] += a*temp[m];
            fshift[    slj][m] += b*temp[m];
        }
    }

    if (VirCorr)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i]*fv[j] + xix[i]*temp[j];
            }
        }
    }

    /* TOTAL: 77 flops */
}


static void spread_vsite4FDN(t_iatom ia[], real a, real b, real c,
                             rvec x[], rvec f[], rvec fshift[],
                             gmx_bool VirCorr, matrix dxdf,
                             t_pbc *pbc, t_graph *g)
{
    rvec xvi, xij, xik, xil, ra, rb, rja, rjb, rab, rm, rt;
    rvec fv, fj, fk, fl;
    real invrm, denom;
    real cfx, cfy, cfz;
    ivec di;
    int  av, ai, aj, ak, al;
    int  svi, sij, sik, sil;

    /* DEBUG: check atom indices */
    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    al = ia[5];

    copy_rvec(f[av], fv);

    sij = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    sik = pbc_rvec_sub(pbc, x[ak], x[ai], xik);
    sil = pbc_rvec_sub(pbc, x[al], x[ai], xil);
    /* 9 flops */

    ra[XX] = a*xik[XX];
    ra[YY] = a*xik[YY];
    ra[ZZ] = a*xik[ZZ];

    rb[XX] = b*xil[XX];
    rb[YY] = b*xil[YY];
    rb[ZZ] = b*xil[ZZ];

    /* 6 flops */

    rvec_sub(ra, xij, rja);
    rvec_sub(rb, xij, rjb);
    rvec_sub(rb, ra, rab);
    /* 9 flops */

    cprod(rja, rjb, rm);
    /* 9 flops */

    invrm = gmx_invsqrt(norm2(rm));
    denom = invrm*invrm;
    /* 5+5+2 flops */

    cfx = c*invrm*fv[XX];
    cfy = c*invrm*fv[YY];
    cfz = c*invrm*fv[ZZ];
    /* 6 Flops */

    cprod(rm, rab, rt);
    /* 9 flops */

    rt[XX] *= denom;
    rt[YY] *= denom;
    rt[ZZ] *= denom;
    /* 3flops */

    fj[XX] = (        -rm[XX]*rt[XX]) * cfx + ( rab[ZZ]-rm[YY]*rt[XX]) * cfy + (-rab[YY]-rm[ZZ]*rt[XX]) * cfz;
    fj[YY] = (-rab[ZZ]-rm[XX]*rt[YY]) * cfx + (        -rm[YY]*rt[YY]) * cfy + ( rab[XX]-rm[ZZ]*rt[YY]) * cfz;
    fj[ZZ] = ( rab[YY]-rm[XX]*rt[ZZ]) * cfx + (-rab[XX]-rm[YY]*rt[ZZ]) * cfy + (        -rm[ZZ]*rt[ZZ]) * cfz;
    /* 30 flops */

    cprod(rjb, rm, rt);
    /* 9 flops */

    rt[XX] *= denom*a;
    rt[YY] *= denom*a;
    rt[ZZ] *= denom*a;
    /* 3flops */

    fk[XX] = (          -rm[XX]*rt[XX]) * cfx + (-a*rjb[ZZ]-rm[YY]*rt[XX]) * cfy + ( a*rjb[YY]-rm[ZZ]*rt[XX]) * cfz;
    fk[YY] = ( a*rjb[ZZ]-rm[XX]*rt[YY]) * cfx + (          -rm[YY]*rt[YY]) * cfy + (-a*rjb[XX]-rm[ZZ]*rt[YY]) * cfz;
    fk[ZZ] = (-a*rjb[YY]-rm[XX]*rt[ZZ]) * cfx + ( a*rjb[XX]-rm[YY]*rt[ZZ]) * cfy + (          -rm[ZZ]*rt[ZZ]) * cfz;
    /* 36 flops */

    cprod(rm, rja, rt);
    /* 9 flops */

    rt[XX] *= denom*b;
    rt[YY] *= denom*b;
    rt[ZZ] *= denom*b;
    /* 3flops */

    fl[XX] = (          -rm[XX]*rt[XX]) * cfx + ( b*rja[ZZ]-rm[YY]*rt[XX]) * cfy + (-b*rja[YY]-rm[ZZ]*rt[XX]) * cfz;
    fl[YY] = (-b*rja[ZZ]-rm[XX]*rt[YY]) * cfx + (          -rm[YY]*rt[YY]) * cfy + ( b*rja[XX]-rm[ZZ]*rt[YY]) * cfz;
    fl[ZZ] = ( b*rja[YY]-rm[XX]*rt[ZZ]) * cfx + (-b*rja[XX]-rm[YY]*rt[ZZ]) * cfy + (          -rm[ZZ]*rt[ZZ]) * cfz;
    /* 36 flops */

    f[ai][XX] += fv[XX] - fj[XX] - fk[XX] - fl[XX];
    f[ai][YY] += fv[YY] - fj[YY] - fk[YY] - fl[YY];
    f[ai][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ] - fl[ZZ];
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    rvec_inc(f[al], fl);
    /* 21 flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, av), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sij = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, ai), di);
        sik = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, al), SHIFT_IVEC(g, ai), di);
        sil = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sij != CENTRAL || sik != CENTRAL || sil != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - fj[XX] - fk[XX] - fl[XX];
        fshift[CENTRAL][YY] += fv[YY] - fj[YY] - fk[YY] - fl[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ] - fl[ZZ];
        rvec_inc(fshift[sij], fj);
        rvec_inc(fshift[sik], fk);
        rvec_inc(fshift[sil], fl);
    }

    if (VirCorr)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i]*fv[j] + xij[i]*fj[j] + xik[i]*fk[j] + xil[i]*fl[j];
            }
        }
    }

    /* Total: 207 flops (Yuck!) */
}


static int spread_vsiten(t_iatom ia[], t_iparams ip[],
                         rvec x[], rvec f[], rvec fshift[],
                         t_pbc *pbc, t_graph *g)
{
    rvec xv, dx, fi;
    int  n3, av, i, ai;
    real a;
    ivec di;
    int  siv;

    n3 = 3*ip[ia[0]].vsiten.n;
    av = ia[1];
    copy_rvec(x[av], xv);

    for (i = 0; i < n3; i += 3)
    {
        ai = ia[i+2];
        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, av), di);
            siv = IVEC2IS(di);
        }
        else if (pbc)
        {
            siv = pbc_dx_aiuc(pbc, x[ai], xv, dx);
        }
        else
        {
            siv = CENTRAL;
        }
        a = ip[ia[i]].vsiten.a;
        svmul(a, f[av], fi);
        rvec_inc(f[ai], fi);
        if (fshift && siv != CENTRAL)
        {
            rvec_inc(fshift[siv], fi);
            rvec_dec(fshift[CENTRAL], fi);
        }
        /* 6 Flops */
    }

    return n3;
}


static int vsite_count(const t_ilist *ilist, int ftype)
{
    if (ftype == F_VSITEN)
    {
        return ilist[ftype].nr/3;
    }
    else
    {
        return ilist[ftype].nr/(1 + interaction_function[ftype].nratoms);
    }
}

static void spread_vsite_f_thread(gmx_vsite_t *vsite,
                                  rvec x[], rvec f[], rvec *fshift,
                                  gmx_bool VirCorr, matrix dxdf,
                                  t_iparams ip[], t_ilist ilist[],
                                  t_graph *g, t_pbc *pbc_null)
{
    gmx_bool   bPBCAll;
    real       a1, b1, c1;
    int        i, inc, nra, nr, tp, ftype;
    t_iatom   *ia;
    t_pbc     *pbc_null2;
    int       *vsite_pbc;

    if (VirCorr)
    {
        clear_mat(dxdf);
    }

    bPBCAll = (pbc_null != NULL && !vsite->bHaveChargeGroups);

    /* this loop goes backwards to be able to build *
     * higher type vsites from lower types         */
    pbc_null2 = NULL;
    vsite_pbc = NULL;
    for (ftype = F_NRE-1; (ftype >= 0); ftype--)
    {
        if ((interaction_function[ftype].flags & IF_VSITE) &&
            ilist[ftype].nr > 0)
        {
            nra    = interaction_function[ftype].nratoms;
            inc    = 1 + nra;
            nr     = ilist[ftype].nr;
            ia     = ilist[ftype].iatoms;

            if (bPBCAll)
            {
                pbc_null2 = pbc_null;
            }
            else if (pbc_null != NULL)
            {
                vsite_pbc = vsite->vsite_pbc_loc[ftype-F_VSITE2];
            }

            for (i = 0; i < nr; )
            {
                if (vsite_pbc != NULL)
                {
                    if (vsite_pbc[i/(1+nra)] > -2)
                    {
                        pbc_null2 = pbc_null;
                    }
                    else
                    {
                        pbc_null2 = NULL;
                    }
                }

                tp   = ia[0];

                /* Constants for constructing */
                a1   = ip[tp].vsite.a;
                /* Construct the vsite depending on type */
                switch (ftype)
                {
                    case F_VSITE2:
                        spread_vsite2(ia, a1, x, f, fshift, pbc_null2, g);
                        break;
                    case F_VSITE3:
                        b1 = ip[tp].vsite.b;
                        spread_vsite3(ia, a1, b1, x, f, fshift, pbc_null2, g);
                        break;
                    case F_VSITE3FD:
                        b1 = ip[tp].vsite.b;
                        spread_vsite3FD(ia, a1, b1, x, f, fshift, VirCorr, dxdf, pbc_null2, g);
                        break;
                    case F_VSITE3FAD:
                        b1 = ip[tp].vsite.b;
                        spread_vsite3FAD(ia, a1, b1, x, f, fshift, VirCorr, dxdf, pbc_null2, g);
                        break;
                    case F_VSITE3OUT:
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        spread_vsite3OUT(ia, a1, b1, c1, x, f, fshift, VirCorr, dxdf, pbc_null2, g);
                        break;
                    case F_VSITE4FD:
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        spread_vsite4FD(ia, a1, b1, c1, x, f, fshift, VirCorr, dxdf, pbc_null2, g);
                        break;
                    case F_VSITE4FDN:
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        spread_vsite4FDN(ia, a1, b1, c1, x, f, fshift, VirCorr, dxdf, pbc_null2, g);
                        break;
                    case F_VSITEN:
                        inc = spread_vsiten(ia, ip, x, f, fshift, pbc_null2, g);
                        break;
                    default:
                        gmx_fatal(FARGS, "No such vsite type %d in %s, line %d",
                                  ftype, __FILE__, __LINE__);
                }
                clear_rvec(f[ia[1]]);

                /* Increment loop variables */
                i  += inc;
                ia += inc;
            }
        }
    }
}

void spread_vsite_f(gmx_vsite_t *vsite,
                    rvec x[], rvec f[], rvec *fshift,
                    gmx_bool VirCorr, matrix vir,
                    t_nrnb *nrnb, t_idef *idef,
                    int ePBC, gmx_bool bMolPBC, t_graph *g, matrix box,
                    t_commrec *cr)
{
    t_pbc pbc, *pbc_null;
    int   th;

    /* We only need to do pbc when we have inter-cg vsites */
    if ((DOMAINDECOMP(cr) || bMolPBC) && vsite->n_intercg_vsite)
    {
        /* This is wasting some CPU time as we now do this multiple times
         * per MD step. But how often do we have vsites with full pbc?
         */
        pbc_null = set_pbc_dd(&pbc, ePBC, cr->dd, FALSE, box);
    }
    else
    {
        pbc_null = NULL;
    }

    if (DOMAINDECOMP(cr))
    {
        dd_clear_f_vsites(cr->dd, f);
    }

    if (vsite->nthreads == 1)
    {
        spread_vsite_f_thread(vsite,
                              x, f, fshift,
                              VirCorr, vsite->tdata[0].dxdf,
                              idef->iparams, idef->il,
                              g, pbc_null);
    }
    else
    {
        /* First spread the vsites that might depend on other vsites */
        spread_vsite_f_thread(vsite,
                              x, f, fshift,
                              VirCorr, vsite->tdata[vsite->nthreads].dxdf,
                              idef->iparams,
                              vsite->tdata[vsite->nthreads].ilist,
                              g, pbc_null);

#pragma omp parallel num_threads(vsite->nthreads)
        {
            int   thread;
            rvec *fshift_t;

            thread = gmx_omp_get_thread_num();

            if (thread == 0 || fshift == NULL)
            {
                fshift_t = fshift;
            }
            else
            {
                int i;

                fshift_t = vsite->tdata[thread].fshift;

                for (i = 0; i < SHIFTS; i++)
                {
                    clear_rvec(fshift_t[i]);
                }
            }

            spread_vsite_f_thread(vsite,
                                  x, f, fshift_t,
                                  VirCorr, vsite->tdata[thread].dxdf,
                                  idef->iparams,
                                  vsite->tdata[thread].ilist,
                                  g, pbc_null);
        }

        if (fshift != NULL)
        {
            int i;

            for (th = 1; th < vsite->nthreads; th++)
            {
                for (i = 0; i < SHIFTS; i++)
                {
                    rvec_inc(fshift[i], vsite->tdata[th].fshift[i]);
                }
            }
        }
    }

    if (VirCorr)
    {
        int i, j;

        for (th = 0; th < (vsite->nthreads == 1 ? 1 : vsite->nthreads+1); th++)
        {
            for (i = 0; i < DIM; i++)
            {
                for (j = 0; j < DIM; j++)
                {
                    vir[i][j] += -0.5*vsite->tdata[th].dxdf[i][j];
                }
            }
        }
    }

    if (DOMAINDECOMP(cr))
    {
        dd_move_f_vsites(cr->dd, f, fshift);
    }

    inc_nrnb(nrnb, eNR_VSITE2,   vsite_count(idef->il, F_VSITE2));
    inc_nrnb(nrnb, eNR_VSITE3,   vsite_count(idef->il, F_VSITE3));
    inc_nrnb(nrnb, eNR_VSITE3FD, vsite_count(idef->il, F_VSITE3FD));
    inc_nrnb(nrnb, eNR_VSITE3FAD, vsite_count(idef->il, F_VSITE3FAD));
    inc_nrnb(nrnb, eNR_VSITE3OUT, vsite_count(idef->il, F_VSITE3OUT));
    inc_nrnb(nrnb, eNR_VSITE4FD, vsite_count(idef->il, F_VSITE4FD));
    inc_nrnb(nrnb, eNR_VSITE4FDN, vsite_count(idef->il, F_VSITE4FDN));
    inc_nrnb(nrnb, eNR_VSITEN,   vsite_count(idef->il, F_VSITEN));
}

static int *atom2cg(t_block *cgs)
{
    int *a2cg, cg, i;

    snew(a2cg, cgs->index[cgs->nr]);
    for (cg = 0; cg < cgs->nr; cg++)
    {
        for (i = cgs->index[cg]; i < cgs->index[cg+1]; i++)
        {
            a2cg[i] = cg;
        }
    }

    return a2cg;
}

static int count_intercg_vsite(gmx_mtop_t *mtop,
                               gmx_bool   *bHaveChargeGroups)
{
    int             mb, ftype, nral, i, cg, a;
    gmx_molblock_t *molb;
    gmx_moltype_t  *molt;
    int            *a2cg;
    t_ilist        *il;
    t_iatom        *ia;
    int             n_intercg_vsite;

    *bHaveChargeGroups = FALSE;

    n_intercg_vsite = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];

        if (molt->cgs.nr < molt->atoms.nr)
        {
            *bHaveChargeGroups = TRUE;
        }

        a2cg = atom2cg(&molt->cgs);
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (interaction_function[ftype].flags & IF_VSITE)
            {
                nral = NRAL(ftype);
                il   = &molt->ilist[ftype];
                ia   = il->iatoms;
                for (i = 0; i < il->nr; i += 1+nral)
                {
                    cg = a2cg[ia[1+i]];
                    for (a = 1; a < nral; a++)
                    {
                        if (a2cg[ia[1+a]] != cg)
                        {
                            n_intercg_vsite += molb->nmol;
                            break;
                        }
                    }
                }
            }
        }
        sfree(a2cg);
    }

    return n_intercg_vsite;
}

static int **get_vsite_pbc(t_iparams *iparams, t_ilist *ilist,
                           t_atom *atom, t_mdatoms *md,
                           t_block *cgs, int *a2cg)
{
    int      ftype, nral, i, j, vsi, vsite, cg_v, a, nc3 = 0;
    t_ilist *il;
    t_iatom *ia;
    int    **vsite_pbc, *vsite_pbc_f;
    char    *pbc_set;
    gmx_bool bViteOnlyCG_and_FirstAtom;

    /* Make an array that tells if the pbc of an atom is set */
    snew(pbc_set, cgs->index[cgs->nr]);
    /* PBC is set for all non vsites */
    for (a = 0; a < cgs->index[cgs->nr]; a++)
    {
        if ((atom && atom[a].ptype != eptVSite) ||
            (md   && md->ptype[a]  != eptVSite))
        {
            pbc_set[a] = 1;
        }
    }

    snew(vsite_pbc, F_VSITEN-F_VSITE2+1);

    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nral = NRAL(ftype);
            il   = &ilist[ftype];
            ia   = il->iatoms;

            snew(vsite_pbc[ftype-F_VSITE2], il->nr/(1+nral));
            vsite_pbc_f = vsite_pbc[ftype-F_VSITE2];

            i = 0;
            while (i < il->nr)
            {
                vsi   = i/(1+nral);
                vsite = ia[i+1];
                cg_v  = a2cg[vsite];
                /* A value of -2 signals that this vsite and its contructing
                 * atoms are all within the same cg, so no pbc is required.
                 */
                vsite_pbc_f[vsi] = -2;
                /* Check if constructing atoms are outside the vsite's cg */
                nc3 = 0;
                if (ftype == F_VSITEN)
                {
                    nc3 = 3*iparams[ia[i]].vsiten.n;
                    for (j = 0; j < nc3; j += 3)
                    {
                        if (a2cg[ia[i+j+2]] != cg_v)
                        {
                            vsite_pbc_f[vsi] = -1;
                        }
                    }
                }
                else
                {
                    for (a = 1; a < nral; a++)
                    {
                        if (a2cg[ia[i+1+a]] != cg_v)
                        {
                            vsite_pbc_f[vsi] = -1;
                        }
                    }
                }
                if (vsite_pbc_f[vsi] == -1)
                {
                    /* Check if this is the first processed atom of a vsite only cg */
                    bViteOnlyCG_and_FirstAtom = TRUE;
                    for (a = cgs->index[cg_v]; a < cgs->index[cg_v+1]; a++)
                    {
                        /* Non-vsites already have pbc set, so simply check for pbc_set */
                        if (pbc_set[a])
                        {
                            bViteOnlyCG_and_FirstAtom = FALSE;
                            break;
                        }
                    }
                    if (bViteOnlyCG_and_FirstAtom)
                    {
                        /* First processed atom of a vsite only charge group.
                         * The pbc of the input coordinates to construct_vsites
                         * should be preserved.
                         */
                        vsite_pbc_f[vsi] = vsite;
                    }
                    else if (cg_v != a2cg[ia[1+i+1]])
                    {
                        /* This vsite has a different charge group index
                         * than it's first constructing atom
                         * and the charge group has more than one atom,
                         * search for the first normal particle
                         * or vsite that already had its pbc defined.
                         * If nothing is found, use full pbc for this vsite.
                         */
                        for (a = cgs->index[cg_v]; a < cgs->index[cg_v+1]; a++)
                        {
                            if (a != vsite && pbc_set[a])
                            {
                                vsite_pbc_f[vsi] = a;
                                if (gmx_debug_at)
                                {
                                    fprintf(debug, "vsite %d match pbc with atom %d\n",
                                            vsite+1, a+1);
                                }
                                break;
                            }
                        }
                        if (gmx_debug_at)
                        {
                            fprintf(debug, "vsite atom %d  cg %d - %d pbc atom %d\n",
                                    vsite+1, cgs->index[cg_v]+1, cgs->index[cg_v+1],
                                    vsite_pbc_f[vsi]+1);
                        }
                    }
                }
                if (ftype == F_VSITEN)
                {
                    /* The other entries in vsite_pbc_f are not used for center vsites */
                    i += nc3;
                }
                else
                {
                    i += 1+nral;
                }

                /* This vsite now has its pbc defined */
                pbc_set[vsite] = 1;
            }
        }
    }

    sfree(pbc_set);

    return vsite_pbc;
}


gmx_vsite_t *init_vsite(gmx_mtop_t *mtop, t_commrec *cr,
                        gmx_bool bSerial_NoPBC)
{
    int            nvsite, i;
    int           *a2cg;
    gmx_vsite_t   *vsite;
    int            mt;
    gmx_moltype_t *molt;

    /* check if there are vsites */
    nvsite = 0;
    for (i = 0; i < F_NRE; i++)
    {
        if (interaction_function[i].flags & IF_VSITE)
        {
            nvsite += gmx_mtop_ftype_count(mtop, i);
        }
    }

    if (nvsite == 0)
    {
        return NULL;
    }

    snew(vsite, 1);

    vsite->n_intercg_vsite = count_intercg_vsite(mtop,
                                                 &vsite->bHaveChargeGroups);

    /* If we don't have charge groups, the vsite follows its own pbc */
    if (!bSerial_NoPBC &&
        vsite->bHaveChargeGroups &&
        vsite->n_intercg_vsite > 0 && DOMAINDECOMP(cr))
    {
        vsite->nvsite_pbc_molt = mtop->nmoltype;
        snew(vsite->vsite_pbc_molt, vsite->nvsite_pbc_molt);
        for (mt = 0; mt < mtop->nmoltype; mt++)
        {
            molt = &mtop->moltype[mt];
            /* Make an atom to charge group index */
            a2cg = atom2cg(&molt->cgs);
            vsite->vsite_pbc_molt[mt] = get_vsite_pbc(mtop->ffparams.iparams,
                                                      molt->ilist,
                                                      molt->atoms.atom, NULL,
                                                      &molt->cgs, a2cg);
            sfree(a2cg);
        }

        snew(vsite->vsite_pbc_loc_nalloc, F_VSITEN-F_VSITE2+1);
        snew(vsite->vsite_pbc_loc, F_VSITEN-F_VSITE2+1);
    }

    if (bSerial_NoPBC)
    {
        vsite->nthreads = 1;
    }
    else
    {
        vsite->nthreads = gmx_omp_nthreads_get(emntVSITE);
    }
    if (!bSerial_NoPBC)
    {
        /* We need one extra thread data structure for the overlap vsites */
        snew(vsite->tdata, vsite->nthreads+1);
    }

    vsite->th_ind        = NULL;
    vsite->th_ind_nalloc = 0;

    return vsite;
}

static void prepare_vsite_thread(const t_ilist      *ilist,
                                 gmx_vsite_thread_t *vsite_th)
{
    int ftype;

    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            if (ilist[ftype].nr > vsite_th->ilist[ftype].nalloc)
            {
                vsite_th->ilist[ftype].nalloc = over_alloc_large(ilist[ftype].nr);
                srenew(vsite_th->ilist[ftype].iatoms, vsite_th->ilist[ftype].nalloc);
            }

            vsite_th->ilist[ftype].nr = 0;
        }
    }
}

void split_vsites_over_threads(const t_ilist   *ilist,
                               const t_iparams *ip,
                               const t_mdatoms *mdatoms,
                               gmx_bool         bLimitRange,
                               gmx_vsite_t     *vsite)
{
    int      th;
    int      vsite_atom_range, natperthread;
    int     *th_ind;
    int      ftype;
    t_iatom *iat;
    t_ilist *il_th;
    int      nral1, inc, i, j;

    if (vsite->nthreads == 1)
    {
        /* Nothing to do */
        return;
    }

#pragma omp parallel for num_threads(vsite->nthreads) schedule(static)
    for (th = 0; th < vsite->nthreads; th++)
    {
        prepare_vsite_thread(ilist, &vsite->tdata[th]);
    }
    /* Master threads does the (potential) overlap vsites */
    prepare_vsite_thread(ilist, &vsite->tdata[vsite->nthreads]);

    /* The current way of distributing the vsites over threads in primitive.
     * We divide the atom range 0 - natoms_in_vsite uniformly over threads,
     * without taking into account how the vsites are distributed.
     * Without domain decomposition we bLimitRange=TRUE and we at least
     * tighten the upper bound of the range (useful for common systems
     * such as a vsite-protein in 3-site water).
     */
    if (bLimitRange)
    {
        vsite_atom_range = -1;
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (interaction_function[ftype].flags & IF_VSITE)
            {
                if (ftype != F_VSITEN)
                {
                    nral1 = 1 + NRAL(ftype);
                    iat   = ilist[ftype].iatoms;
                    for (i = 0; i < ilist[ftype].nr; i += nral1)
                    {
                        for (j = i + 1; j < i + nral1; j++)
                        {
                            vsite_atom_range = std::max(vsite_atom_range, iat[j]);
                        }
                    }
                }
                else
                {
                    int vs_ind_end;

                    iat = ilist[ftype].iatoms;

                    i = 0;
                    while (i < ilist[ftype].nr)
                    {
                        /* The 3 below is from 1+NRAL(ftype)=3 */
                        vs_ind_end = i + ip[iat[i]].vsiten.n*3;

                        vsite_atom_range = std::max(vsite_atom_range, iat[i+1]);
                        while (i < vs_ind_end)
                        {
                            vsite_atom_range = std::max(vsite_atom_range, iat[i+2]);
                            i               += 3;
                        }
                    }
                }
            }
        }
        vsite_atom_range++;
    }
    else
    {
        vsite_atom_range = mdatoms->homenr;
    }
    natperthread = (vsite_atom_range + vsite->nthreads - 1)/vsite->nthreads;

    if (debug)
    {
        fprintf(debug, "virtual site thread dist: natoms %d, range %d, natperthread %d\n", mdatoms->nr, vsite_atom_range, natperthread);
    }

    /* To simplify the vsite assignment, we make an index which tells us
     * to which thread particles, both non-vsites and vsites, are assigned.
     */
    if (mdatoms->nr > vsite->th_ind_nalloc)
    {
        vsite->th_ind_nalloc = over_alloc_large(mdatoms->nr);
        srenew(vsite->th_ind, vsite->th_ind_nalloc);
    }
    th_ind = vsite->th_ind;
    th     = 0;
    for (i = 0; i < mdatoms->nr; i++)
    {
        if (mdatoms->ptype[i] == eptVSite)
        {
            /* vsites are not assigned to a thread yet */
            th_ind[i] = -1;
        }
        else
        {
            /* assign non-vsite particles to thread th */
            th_ind[i] = th;
        }
        if (i == (th + 1)*natperthread && th < vsite->nthreads)
        {
            th++;
        }
    }

    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nral1 = 1 + NRAL(ftype);
            inc   = nral1;
            iat   = ilist[ftype].iatoms;
            for (i = 0; i < ilist[ftype].nr; )
            {
                th = iat[1+i]/natperthread;
                /* We would like to assign this vsite the thread th,
                 * but it might depend on atoms outside the atom range of th
                 * or on another vsite not assigned to thread th.
                 */
                if (ftype != F_VSITEN)
                {
                    for (j = i + 2; j < i + nral1; j++)
                    {
                        if (th_ind[iat[j]] != th)
                        {
                            /* Some constructing atoms are not assigned to
                             * thread th, move this vsite to a separate batch.
                             */
                            th = vsite->nthreads;
                        }
                    }
                }
                else
                {
                    /* The 3 below is from 1+NRAL(ftype)=3 */
                    inc = ip[iat[i]].vsiten.n*3;
                    for (j = i + 2; j < i + inc; j += 3)
                    {
                        if (th_ind[iat[j]] != th)
                        {
                            th = vsite->nthreads;
                        }
                    }
                }
                /* Copy this vsite to the thread data struct of thread th */
                il_th = &vsite->tdata[th].ilist[ftype];
                for (j = i; j < i + inc; j++)
                {
                    il_th->iatoms[il_th->nr++] = iat[j];
                }
                /* Update this vsite's thread index entry */
                th_ind[iat[1+i]] = th;

                i += inc;
            }
        }
    }

    if (debug)
    {
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if ((interaction_function[ftype].flags & IF_VSITE) &&
                ilist[ftype].nr > 0)
            {
                fprintf(debug, "%-20s thread dist:",
                        interaction_function[ftype].longname);
                for (th = 0; th < vsite->nthreads+1; th++)
                {
                    fprintf(debug, " %4d", vsite->tdata[th].ilist[ftype].nr);
                }
                fprintf(debug, "\n");
            }
        }
    }
}

void set_vsite_top(gmx_vsite_t *vsite, gmx_localtop_t *top, t_mdatoms *md,
                   t_commrec *cr)
{
    int *a2cg;

    if (vsite->n_intercg_vsite > 0)
    {
        if (vsite->bHaveChargeGroups)
        {
            /* Make an atom to charge group index */
            a2cg                 = atom2cg(&top->cgs);
            vsite->vsite_pbc_loc = get_vsite_pbc(top->idef.iparams,
                                                 top->idef.il, NULL, md,
                                                 &top->cgs, a2cg);
            sfree(a2cg);
        }

    }

    if (vsite->nthreads > 1)
    {
        if (vsite->bHaveChargeGroups)
        {
            gmx_fatal(FARGS, "The combination of threading, virtual sites and charge groups is not implemented");
        }

        split_vsites_over_threads(top->idef.il, top->idef.iparams,
                                  md, !DOMAINDECOMP(cr), vsite);
    }
}
