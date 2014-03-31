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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "gromacs/utility/smalloc.h"
#include "pbc.h"
#include "txtdump.h"
#include "vec.h"
#include "nrnb.h"
#include "constr.h"

typedef struct gmx_shakedata
{
    rvec *rij;
    real *M2;
    real *tt;
    real *dist2;
    int   nalloc;
    /* SOR stuff */
    real  delta;
    real  omega;
    real  gamma;
} t_gmx_shakedata;

gmx_shakedata_t shake_init()
{
    gmx_shakedata_t d;

    snew(d, 1);

    d->nalloc = 0;
    d->rij    = NULL;
    d->M2     = NULL;
    d->tt     = NULL;
    d->dist2  = NULL;

    /* SOR initialization */
    d->delta = 0.1;
    d->omega = 1.0;
    d->gamma = 1000000;

    return d;
}

static void pv(FILE *log, char *s, rvec x)
{
    int m;

    fprintf(log, "%5s:", s);
    for (m = 0; (m < DIM); m++)
    {
        fprintf(log, "  %10.3f", x[m]);
    }
    fprintf(log, "\n");
    fflush(log);
}

void cshake(atom_id iatom[], int ncon, int *nnit, int maxnit,
            real dist2[], real xp[], real rij[], real m2[], real omega,
            real invmass[], real tt[], real lagr[], int *nerror)
{
    /*
     *     r.c. van schaik and w.f. van gunsteren
     *     eth zuerich
     *     june 1992
     *     Adapted for use with Gromacs by David van der Spoel november 92 and later.
     */
    /* default should be increased! MRS 8/4/2009 */
    const   real mytol = 1e-10;

    int          ll, i, j, i3, j3, l3;
    int          ix, iy, iz, jx, jy, jz;
    real         toler, rpij2, rrpr, tx, ty, tz, diff, acor, im, jm;
    real         xh, yh, zh, rijx, rijy, rijz;
    real         tix, tiy, tiz;
    real         tjx, tjy, tjz;
    int          nit, error, nconv;
    real         iconvf;

    error = 0;
    nconv = 1;
    for (nit = 0; (nit < maxnit) && (nconv != 0) && (error == 0); nit++)
    {
        nconv = 0;
        for (ll = 0; (ll < ncon) && (error == 0); ll++)
        {
            l3    = 3*ll;
            rijx  = rij[l3+XX];
            rijy  = rij[l3+YY];
            rijz  = rij[l3+ZZ];
            i     = iatom[l3+1];
            j     = iatom[l3+2];
            i3    = 3*i;
            j3    = 3*j;
            ix    = i3+XX;
            iy    = i3+YY;
            iz    = i3+ZZ;
            jx    = j3+XX;
            jy    = j3+YY;
            jz    = j3+ZZ;

            tx      = xp[ix]-xp[jx];
            ty      = xp[iy]-xp[jy];
            tz      = xp[iz]-xp[jz];
            rpij2   = tx*tx+ty*ty+tz*tz;
            toler   = dist2[ll];
            diff    = toler-rpij2;

            /* iconvf is less than 1 when the error is smaller than a bound */
            /* But if tt is too big, then it will result in looping in iconv */

            iconvf = fabs(diff)*tt[ll];

            if (iconvf > 1)
            {
                nconv   = iconvf;
                rrpr    = rijx*tx+rijy*ty+rijz*tz;

                if (rrpr < toler*mytol)
                {
                    error = ll+1;
                }
                else
                {
                    acor      = omega*diff*m2[ll]/rrpr;
                    lagr[ll] += acor;
                    xh        = rijx*acor;
                    yh        = rijy*acor;
                    zh        = rijz*acor;
                    im        = invmass[i];
                    jm        = invmass[j];
                    xp[ix]   += xh*im;
                    xp[iy]   += yh*im;
                    xp[iz]   += zh*im;
                    xp[jx]   -= xh*jm;
                    xp[jy]   -= yh*jm;
                    xp[jz]   -= zh*jm;
                }
            }
        }
    }
    *nnit   = nit;
    *nerror = error;
}

int vec_shakef(FILE *fplog, gmx_shakedata_t shaked,
               real invmass[], int ncon,
               t_iparams ip[], t_iatom *iatom,
               real tol, rvec x[], rvec prime[], real omega,
               gmx_bool bFEP, real lambda, real lagr[],
               real invdt, rvec *v,
               gmx_bool bCalcVir, tensor vir_r_m_dr, int econq,
               t_vetavars *vetavar)
{
    rvec    *rij;
    real    *M2, *tt, *dist2;
    int      maxnit = 1000;
    int      nit    = 0, ll, i, j, type;
    t_iatom *ia;
    real     L1, tol2, toler;
    real     mm    = 0., tmp;
    int      error = 0;
    real     g, vscale, rscale, rvscale;

    if (ncon > shaked->nalloc)
    {
        shaked->nalloc = over_alloc_dd(ncon);
        srenew(shaked->rij, shaked->nalloc);
        srenew(shaked->M2, shaked->nalloc);
        srenew(shaked->tt, shaked->nalloc);
        srenew(shaked->dist2, shaked->nalloc);
    }
    rij   = shaked->rij;
    M2    = shaked->M2;
    tt    = shaked->tt;
    dist2 = shaked->dist2;

    L1   = 1.0-lambda;
    tol2 = 2.0*tol;
    ia   = iatom;
    for (ll = 0; (ll < ncon); ll++, ia += 3)
    {
        type  = ia[0];
        i     = ia[1];
        j     = ia[2];

        mm          = 2*(invmass[i]+invmass[j]);
        rij[ll][XX] = x[i][XX]-x[j][XX];
        rij[ll][YY] = x[i][YY]-x[j][YY];
        rij[ll][ZZ] = x[i][ZZ]-x[j][ZZ];
        M2[ll]      = 1.0/mm;
        if (bFEP)
        {
            toler = sqr(L1*ip[type].constr.dA + lambda*ip[type].constr.dB);
        }
        else
        {
            toler = sqr(ip[type].constr.dA);
        }
        dist2[ll] = toler;
        tt[ll]    = 1.0/(toler*tol2);
    }

    switch (econq)
    {
        case econqCoord:
            cshake(iatom, ncon, &nit, maxnit, dist2, prime[0], rij[0], M2, omega, invmass, tt, lagr, &error);
            break;
        case econqVeloc:
            crattle(iatom, ncon, &nit, maxnit, dist2, prime[0], rij[0], M2, omega, invmass, tt, lagr, &error, invdt, vetavar);
            break;
    }

    if (nit >= maxnit)
    {
        if (fplog)
        {
            fprintf(fplog, "Shake did not converge in %d steps\n", maxnit);
        }
        fprintf(stderr, "Shake did not converge in %d steps\n", maxnit);
        nit = 0;
    }
    else if (error != 0)
    {
        if (fplog)
        {
            fprintf(fplog, "Inner product between old and new vector <= 0.0!\n"
                    "constraint #%d atoms %u and %u\n",
                    error-1, iatom[3*(error-1)+1]+1, iatom[3*(error-1)+2]+1);
        }
        fprintf(stderr, "Inner product between old and new vector <= 0.0!\n"
                "constraint #%d atoms %u and %u\n",
                error-1, iatom[3*(error-1)+1]+1, iatom[3*(error-1)+2]+1);
        nit = 0;
    }

    /* Constraint virial and correct the lagrange multipliers for the length */

    ia = iatom;

    for (ll = 0; (ll < ncon); ll++, ia += 3)
    {

        if ((econq == econqCoord) && v != NULL)
        {
            /* Correct the velocities */
            mm = lagr[ll]*invmass[ia[1]]*invdt/vetavar->rscale;
            for (i = 0; i < DIM; i++)
            {
                v[ia[1]][i] += mm*rij[ll][i];
            }
            mm = lagr[ll]*invmass[ia[2]]*invdt/vetavar->rscale;
            for (i = 0; i < DIM; i++)
            {
                v[ia[2]][i] -= mm*rij[ll][i];
            }
            /* 16 flops */
        }

        /* constraint virial */
        if (bCalcVir)
        {
            if (econq == econqCoord)
            {
                mm = lagr[ll]/vetavar->rvscale;
            }
            if (econq == econqVeloc)
            {
                mm = lagr[ll]/(vetavar->vscale*vetavar->vscale_nhc[0]);
            }
            for (i = 0; i < DIM; i++)
            {
                tmp = mm*rij[ll][i];
                for (j = 0; j < DIM; j++)
                {
                    vir_r_m_dr[i][j] -= tmp*rij[ll][j];
                }
            }
            /* 21 flops */
        }

        /* Correct the lagrange multipliers for the length  */
        /* (more details would be useful here . . . )*/

        type  = ia[0];
        if (bFEP)
        {
            toler = L1*ip[type].constr.dA + lambda*ip[type].constr.dB;
        }
        else
        {
            toler     = ip[type].constr.dA;
            lagr[ll] *= toler;
        }
    }

    return nit;
}

static void check_cons(FILE *log, int nc, rvec x[], rvec prime[], rvec v[],
                       t_iparams ip[], t_iatom *iatom,
                       real invmass[], int econq)
{
    t_iatom *ia;
    int      ai, aj;
    int      i;
    real     d, dp;
    rvec     dx, dv;

    fprintf(log,
            "    i     mi      j     mj      before       after   should be\n");
    ia = iatom;
    for (i = 0; (i < nc); i++, ia += 3)
    {
        ai = ia[1];
        aj = ia[2];
        rvec_sub(x[ai], x[aj], dx);
        d = norm(dx);

        switch (econq)
        {
            case econqCoord:
                rvec_sub(prime[ai], prime[aj], dx);
                dp = norm(dx);
                fprintf(log, "%5d  %5.2f  %5d  %5.2f  %10.5f  %10.5f  %10.5f\n",
                        ai+1, 1.0/invmass[ai],
                        aj+1, 1.0/invmass[aj], d, dp, ip[ia[0]].constr.dA);
                break;
            case econqVeloc:
                rvec_sub(v[ai], v[aj], dv);
                d = iprod(dx, dv);
                rvec_sub(prime[ai], prime[aj], dv);
                dp = iprod(dx, dv);
                fprintf(log, "%5d  %5.2f  %5d  %5.2f  %10.5f  %10.5f  %10.5f\n",
                        ai+1, 1.0/invmass[ai],
                        aj+1, 1.0/invmass[aj], d, dp, 0.);
                break;
        }
    }
}

gmx_bool bshakef(FILE *log, gmx_shakedata_t shaked,
                 real invmass[], int nblocks, int sblock[],
                 t_idef *idef, t_inputrec *ir, rvec x_s[], rvec prime[],
                 t_nrnb *nrnb, real *lagr, real lambda, real *dvdlambda,
                 real invdt, rvec *v, gmx_bool bCalcVir, tensor vir_r_m_dr,
                 gmx_bool bDumpOnError, int econq, t_vetavars *vetavar)
{
    t_iatom *iatoms;
    real    *lam, dt_2, dvdl;
    int      i, n0, ncons, blen, type;
    int      tnit = 0, trij = 0;

#ifdef DEBUG
    fprintf(log, "nblocks=%d, sblock[0]=%d\n", nblocks, sblock[0]);
#endif

    ncons = idef->il[F_CONSTR].nr/3;

    for (i = 0; i < ncons; i++)
    {
        lagr[i] = 0;
    }

    iatoms = &(idef->il[F_CONSTR].iatoms[sblock[0]]);
    lam    = lagr;
    for (i = 0; (i < nblocks); )
    {
        blen  = (sblock[i+1]-sblock[i]);
        blen /= 3;
        n0    = vec_shakef(log, shaked, invmass, blen, idef->iparams,
                           iatoms, ir->shake_tol, x_s, prime, shaked->omega,
                           ir->efep != efepNO, lambda, lam, invdt, v, bCalcVir, vir_r_m_dr,
                           econq, vetavar);

#ifdef DEBUGSHAKE
        check_cons(log, blen, x_s, prime, v, idef->iparams, iatoms, invmass, econq);
#endif

        if (n0 == 0)
        {
            if (bDumpOnError && log)
            {
                {
                    check_cons(log, blen, x_s, prime, v, idef->iparams, iatoms, invmass, econq);
                }
            }
            return FALSE;
        }
        tnit   += n0*blen;
        trij   += blen;
        iatoms += 3*blen; /* Increment pointer! */
        lam    += blen;
        i++;
    }
    /* only for position part? */
    if (econq == econqCoord)
    {
        if (ir->efep != efepNO)
        {
            real bondA, bondB;
            dt_2 = 1/sqr(ir->delta_t);
            dvdl = 0;
            for (i = 0; i < ncons; i++)
            {
                type  = idef->il[F_CONSTR].iatoms[3*i];

                /* dh/dl contribution from constraint force is  dh/dr (constraint force) dot dr/dl */
                /* constraint force is -\sum_i lagr_i* d(constraint)/dr, with constrant = r^2-d^2  */
                /* constraint force is -\sum_i lagr_i* 2 r  */
                /* so dh/dl = -\sum_i lagr_i* 2 r * dr/dl */
                /* However, by comparison with lincs and with
                   comparison with a full thermodynamics cycle (see
                   redmine issue #1255), this is off by a factor of
                   two -- the 2r should apparently just be r.  Further
                   investigation should be done at some point to
                   understand why and see if there is something deeper
                   we are missing */

                bondA = idef->iparams[type].constr.dA;
                bondB = idef->iparams[type].constr.dB;
                dvdl += lagr[i] * dt_2 * ((1.0-lambda)*bondA + lambda*bondB) * (bondB-bondA);
            }
            *dvdlambda += dvdl;
        }
    }
#ifdef DEBUG
    fprintf(log, "tnit: %5d  omega: %10.5f\n", tnit, omega);
#endif
    if (ir->bShakeSOR)
    {
        if (tnit > shaked->gamma)
        {
            shaked->delta *= -0.5;
        }
        shaked->omega += shaked->delta;
        shaked->gamma  = tnit;
    }
    inc_nrnb(nrnb, eNR_SHAKE, tnit);
    inc_nrnb(nrnb, eNR_SHAKE_RIJ, trij);
    if (v)
    {
        inc_nrnb(nrnb, eNR_CONSTR_V, trij*2);
    }
    if (bCalcVir)
    {
        inc_nrnb(nrnb, eNR_CONSTR_VIR, trij);
    }

    return TRUE;
}

void crattle(atom_id iatom[], int ncon, int *nnit, int maxnit,
             real dist2[], real vp[], real rij[], real m2[], real omega,
             real invmass[], real tt[], real lagr[], int *nerror, real invdt, t_vetavars *vetavar)
{
    /*
     *     r.c. van schaik and w.f. van gunsteren
     *     eth zuerich
     *     june 1992
     *     Adapted for use with Gromacs by David van der Spoel november 92 and later.
     *     rattle added by M.R. Shirts, April 2004, from code written by Jay Ponder in TINKER
     *     second part of rattle algorithm
     */

    const   real mytol = 1e-10;

    int          ll, i, j, i3, j3, l3, ii;
    int          ix, iy, iz, jx, jy, jz;
    real         toler, rijd, vpijd, vx, vy, vz, diff, acor, xdotd, fac, im, jm, imdt, jmdt;
    real         xh, yh, zh, rijx, rijy, rijz;
    real         tix, tiy, tiz;
    real         tjx, tjy, tjz;
    int          nit, error, nconv;
    real         veta, vscale_nhc, iconvf;

    veta       = vetavar->veta;
    vscale_nhc = vetavar->vscale_nhc[0];  /* for now, just use the first state */

    error = 0;
    nconv = 1;
    for (nit = 0; (nit < maxnit) && (nconv != 0) && (error == 0); nit++)
    {
        nconv = 0;
        for (ll = 0; (ll < ncon) && (error == 0); ll++)
        {
            l3      = 3*ll;
            rijx    = rij[l3+XX];
            rijy    = rij[l3+YY];
            rijz    = rij[l3+ZZ];
            i       = iatom[l3+1];
            j       = iatom[l3+2];
            i3      = 3*i;
            j3      = 3*j;
            ix      = i3+XX;
            iy      = i3+YY;
            iz      = i3+ZZ;
            jx      = j3+XX;
            jy      = j3+YY;
            jz      = j3+ZZ;
            vx      = vp[ix]-vp[jx];
            vy      = vp[iy]-vp[jy];
            vz      = vp[iz]-vp[jz];

            vpijd   = vx*rijx+vy*rijy+vz*rijz;
            toler   = dist2[ll];
            /* this is r(t+dt) \dotproduct \dot{r}(t+dt) */
            xdotd   = vpijd*vscale_nhc + veta*toler;

            /* iconv is zero when the error is smaller than a bound */
            iconvf   = fabs(xdotd)*(tt[ll]/invdt);

            if (iconvf > 1)
            {
                nconv     = iconvf;
                fac       = omega*2.0*m2[ll]/toler;
                acor      = -fac*xdotd;
                lagr[ll] += acor;

                xh        = rijx*acor;
                yh        = rijy*acor;
                zh        = rijz*acor;

                im        = invmass[i]/vscale_nhc;
                jm        = invmass[j]/vscale_nhc;

                vp[ix] += xh*im;
                vp[iy] += yh*im;
                vp[iz] += zh*im;
                vp[jx] -= xh*jm;
                vp[jy] -= yh*jm;
                vp[jz] -= zh*jm;
            }
        }
    }
    *nnit   = nit;
    *nerror = error;
}
